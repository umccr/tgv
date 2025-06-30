use crate::{
    alignment::{Alignment, AlignmentBuilder}, contig::Contig, error::TGVError, helpers::is_url, reference::Reference, region::Region, repository, sequence::Sequence, settings::{BackendType, Settings}, track_service::{TrackService, TrackServiceEnum, UcscDbTrackService}
};
use noodles::bam::io::indexed_reader;
use std::io::{Read, Seek};

use noodles::bam::bai;
use noodles::bam::io::IndexedReader;
use noodles::core::region;
use noodles::core::Position;
use noodles::sam::Header;
use opendal::services::{Gcs, Http, S3};
use opendal::{BlockingOperator, Operator};
use reqwest::Client;
use serde::Deserialize;
use std::path::Path;

pub struct Repository {
    pub alignment_repository: AlignmentRepositoryEnum,
    pub track_service: Option<TrackServiceEnum>,
    pub sequence_service: Option<SequenceService>,
}

impl Repository {
    pub async fn new(settings: &Settings) -> Result<Self, TGVError> {
        let alignment_repository = AlignmentRepositoryEnum::from(settings)?;

        let (track_service, sequence_service): (Option<TrackServiceEnum>, Option<SequenceService>) =
            match settings.reference.as_ref() {
                Some(reference) => {
                    let ts = match settings.backend {
                        //BackendType::Api => TrackServiceEnum::Api(UcscApiTrackService::new()?),
                        BackendType::Db => {
                            TrackServiceEnum::Db(UcscDbTrackService::new(reference).await?)
                        }
                    };
                    let ss = SequenceService::new(reference.clone())?;
                    (Some(ts), Some(ss))
                }
                None => (None, None),
            };

        Ok(Self {
            alignment_repository,
            track_service,
            sequence_service,
        })
    }

    pub fn track_service_checked(&self) -> Result<&TrackServiceEnum, TGVError> {
        match self.track_service {
            Some(ref track_service) => Ok(track_service),
            None => Err(TGVError::StateError(
                "Track service is not initialized".to_string(),
            )),
        }
    }

    pub fn sequence_service_checked(&self) -> Result<&SequenceService, TGVError> {
        match self.sequence_service {
            Some(ref sequence_service) => Ok(sequence_service),
            None => Err(TGVError::StateError(
                "Sequence service is not initialized".to_string(),
            )),
        }
    }

    pub async fn close(&mut self) -> Result<(), TGVError> {
        if let Some(ts) = self.track_service.as_mut() {
            ts.close().await?;
        }
        if let Some(ss) = self.sequence_service.as_mut() {
            ss.close().await?;
        }
        Ok(())
    }

    pub fn has_alignment(&self) -> bool {
        self.alignment_repository.has_alignment()
    }

    pub fn has_track(&self) -> bool {
        self.track_service.is_some()
    }

    pub fn has_sequence(&self) -> bool {
        self.sequence_service.is_some()
    }
}

#[derive(Debug)]
struct RemoteSource {
    op: BlockingOperator,
    key: String,
}

impl RemoteSource {
    fn from(path: &String) -> Result<Self, TGVError> {
        if path.starts_with("s3://") {
            let path = path.strip_prefix("s3://").unwrap();
            let bucket = path.split('/').next().unwrap();
            let key = path.split('/').rev().next().unwrap();
            let root = path
                .strip_prefix(bucket)
                .unwrap()
                .strip_suffix(key)
                .unwrap();

            let builder = S3::default().root(root).bucket(bucket);
            let op = Operator::new(builder).unwrap().finish().blocking();

            Ok(Self {
                op,
                key: key.to_string(),
            })
        } else if path.starts_with("http://") {
            let path = path.strip_prefix("http://").unwrap();
            let endpoint = path.split('/').next().unwrap();
            let key = path.split('/').rev().next().unwrap();
            let root = path
                .strip_prefix(endpoint)
                .unwrap()
                .strip_suffix(key)
                .unwrap();

            let builder = Http::default()
                .endpoint(&["http://".to_string(), endpoint.to_string()].concat())
                .root(root);

            let op = Operator::new(builder).unwrap().finish().blocking();

            Ok(Self {
                op,
                key: key.to_string(),
            })
        } else if path.starts_with("https://") {
            let path = path.strip_prefix("https://").unwrap();
            let endpoint = path.split('/').next().unwrap();
            let key = path.split('/').rev().next().unwrap();
            let root = path
                .strip_prefix(endpoint)
                .unwrap()
                .strip_suffix(key)
                .unwrap();

            let builder = Http::default()
                .endpoint(&["https://".to_string(), endpoint.to_string()].concat())
                .root(root);

            let op = Operator::new(builder).unwrap().finish().blocking();

            Ok(Self {
                op,
                key: key.to_string(),
            })
        } else if path.starts_with("gss://") {
            let path = path.strip_prefix("gss://").unwrap();
            let bucket = path.split('/').next().unwrap();
            let key = path.split('/').rev().next().unwrap();
            let root = path
                .strip_prefix(bucket)
                .unwrap()
                .strip_suffix(key)
                .unwrap();

            let builder = Gcs::default().root(root).bucket(bucket);

            let op = Operator::new(builder).unwrap().finish().blocking();

            Ok(Self {
                op,
                key: key.to_string(),
            })
        } else {
            Err(TGVError::ValueError(format!(
                "Unsupported remote path {}. Only S3, HTTP/HTTPS, and GS are supported.",
                path
            )))
        }
    }

    fn read(self) -> impl Read + Seek {
        let reader = self
            .op
            .reader(&self.key)
            .unwrap()
            .into_std_read(..)
            .unwrap();

        reader
    }
}

pub trait AlignmentRepository {
    fn read_alignment(&self, region: &Region) -> Result<Alignment, TGVError>;

    fn read_header(&self) -> Result<Vec<(String, Option<usize>)>, TGVError>;
}

#[derive(Debug)]
pub struct AlignmentsRepository {
    path: String,
    index_path: Option<String>,
}

impl AlignmentsRepository {
    fn new(path: String, index_path: Option<String>) -> Result<Self, TGVError> {
        if is_url(&path) {
            return Err(TGVError::IOError(format!(
                "{} is a remote path. Use RemoteAlignmentsRepository for remote reads IO",
                path
            )));
        }

        if !Path::new(&path).exists() {
            return Err(TGVError::IOError(format!(
                "Alignments file {} not found",
                path
            )));
        }

        match &index_path {
            Some(index_path) => {
                if !Path::new(index_path).exists() {
                    return Err(TGVError::IOError(format!(
                        "BAM index file {} not found. Only indexed BAM files are supported.",
                        index_path
                    )));
                }
            }
            None => {
                if !Path::new(&format!("{}.bai", path)).exists() {
                    return Err(TGVError::IOError(format!(
                        "BAM index file {}.bai not found. Only indexed BAM files are supported.",
                        path
                    )));
                }
            }
        }

        Ok(Self { path, index_path })
    }
}

impl AlignmentRepository for AlignmentsRepository {
    fn read_alignment(&self, region: &Region) -> Result<Alignment, TGVError> {
        let mut reader = match self.index_path.as_ref() {
            Some(index_path) => {
                let index = bai::fs::read(index_path)?;
                indexed_reader::Builder::default()
                    .set_index(index)
                    .build_from_path(self.path.clone())?
            }
            None => indexed_reader::Builder::default().build_from_path(self.path.clone())?,
        };

        let header = reader.read_header()?;
        let noodles_region = region::Region::new(
            region.contig.name.to_string(),
            Position::new(region.start)
                .ok_or_else(|| TGVError::ValueError("invalid position".to_string()))?
                ..=Position::new(region.end)
                    .ok_or_else(|| TGVError::ValueError("invalid position".to_string()))?,
        );

        let records = reader.query(&header, &noodles_region)?;

        let mut alignment_builder = AlignmentBuilder::new()?;
        for record in records {
            let read = record?;
            alignment_builder.add_read(read)?;
        }

        alignment_builder.region(region)?.build()
    }

    /// Read BAM headers and return contig namesa and lengths.
    /// Note that this function does not interprete the contig name as contg vs chromosome.
    fn read_header(&self) -> Result<Vec<(String, Option<usize>)>, TGVError> {
        let mut reader = match self.index_path.as_ref() {
            Some(index_path) => {
                let index = bai::fs::read(index_path)?;
                indexed_reader::Builder::default()
                    .set_index(index)
                    .build_from_path(self.path.clone())?
            }
            None => indexed_reader::Builder::default().build_from_path(self.path.clone())?,
        };

        let header = reader.read_header()?;
        get_contig_names_and_lengths_from_header(&header)
    }
}

#[derive(Debug)]
pub struct RemoteAlignmentsRepository {
    path: String,
    source: RemoteSource,
}

impl AlignmentRepository for RemoteAlignmentsRepository {
    fn read_alignment(&self, region: &Region) -> Result<Alignment, TGVError> {
        let index = RemoteSource::from(&[&self.path, ".bai"].concat())?.read();
        let source = RemoteSource::from(&self.path)?.read();

        let mut index_reader = bai::io::Reader::new(index);
        let index = index_reader.read_index()?;

        let mut bam = IndexedReader::new(source, index);

        let header = bam.read_header()?;

        let noodles_region = region::Region::new(
            region.contig.name.to_string(),
            Position::new(region.start)
                .ok_or_else(|| TGVError::ValueError("invalid position".to_string()))?
                ..=Position::new(region.end)
                    .ok_or_else(|| TGVError::ValueError("invalid position".to_string()))?,
        );
        let records = bam.query(&header, &noodles_region)?;

        let mut alignment_builder = AlignmentBuilder::new()?;
        for record in records {
            let read = record?;
            alignment_builder.add_read(read)?;
        }

        alignment_builder.region(region)?.build()
        // let query_contig_string = get_query_contig_string(&header, region)?;
        // bam.fetch((
        //     &query_contig_string,
        //     region.start as i32 - 1,
        //     region.end as i32,
        // ))
        // .map_err(|e| TGVError::IOError(e.to_string()))?;

        // let mut alignment_builder = AlignmentBuilder::new()?;
        //
        // for record in bam.records() {
        //     let read = record.map_err(|e| TGVError::IOError(e.to_string()))?;
        //     alignment_builder.add_read(read)?;
        // }
        //
        // alignment_builder.region(region)?.build()
    }

    fn read_header(&self) -> Result<Vec<(String, Option<usize>)>, TGVError> {
        let index = RemoteSource::from(&[&self.path, ".bai"].concat())?.read();
        let source = RemoteSource::from(&self.path)?.read();

        let mut index_reader = bai::io::Reader::new(index);
        let index = index_reader.read_index()?;

        let mut bam = IndexedReader::new(source, index);

        let header = bam.read_header()?;
        get_contig_names_and_lengths_from_header(&header)
    }
}

// fn is_remote_path {
//     IndexedReader::from_url(
//         &Url::parse(path).map_err(|e| TGVError::IOError(e.to_string()))?,
//     )
//     .unwrap();

// struct CRAMRepository {
//     cram_path: String,
// }

fn get_contig_names_and_lengths_from_header(
    header: &Header,
) -> Result<Vec<(String, Option<usize>)>, TGVError> {
    let mut output = Vec::new();

    for (key, reference_sequence) in header.reference_sequences() {
        let length = reference_sequence.length();
        output.push((key.to_string(), Some(length.get())));
    }

    Ok(output)
}

#[derive(Debug)]
pub enum AlignmentRepositoryEnum {
    None,
    OpenDAL(String)
}

impl AlignmentRepositoryEnum {
    pub fn from(settings: &Settings) -> Result<Self, TGVError> {
        if settings.path.is_none() {
            return Ok(AlignmentRepositoryEnum::None);
        }

        let path = settings.path.clone().unwrap();

        if is_url(&path) {
            return Ok(AlignmentRepositoryEnum::OpenDAL(path))
        }

        Ok(AlignmentRepositoryEnum::OpenDAL(path))
    }

    pub fn has_alignment(&self) -> bool {
        match self {
            AlignmentRepositoryEnum::OpenDAL(_path) => true,
            AlignmentRepositoryEnum::None => false,
        }
    }
}

impl AlignmentRepository for AlignmentRepositoryEnum {
    fn read_alignment(&self, region: &Region) -> Result<Alignment, TGVError> {
        match self {
            AlignmentRepositoryEnum::OpenDAL(path) => {
                if is_url(path) {
                    let repo = RemoteAlignmentsRepository {
                        path: path.clone(),
                        source: RemoteSource::from(path)?,
                    };
                    repo.read_alignment(region)
                } else {
                    let repo = AlignmentsRepository::new(path.clone(), None)?;
                    repo.read_alignment(region)
                }
            },
            AlignmentRepositoryEnum::None => Err(TGVError::IOError("No alignment".to_string())),
        }
    }

    fn read_header(&self) -> Result<Vec<(String, Option<usize>)>, TGVError> {
        match self {
            AlignmentRepositoryEnum::OpenDAL(_) => {
                match self {
                    AlignmentRepositoryEnum::OpenDAL(path) => {
                        if is_url(path) {
                            let repo = RemoteAlignmentsRepository {
                                path: path.clone(),
                                source: RemoteSource::from(path)?,
                            };
                            repo.read_header()
                        } else {
                            let repo = AlignmentsRepository::new(path.clone(), None)?;
                            repo.read_header()
                        }
                    },
                    AlignmentRepositoryEnum::None => Err(TGVError::IOError("No alignment".to_string())),
                }
            },
            AlignmentRepositoryEnum::None => Err(TGVError::IOError("No alignment".to_string())),
        }
    }
}

#[derive(Debug, Deserialize)]
struct UcscResponse {
    dna: String,
}

pub struct SequenceService {
    client: Client,
    reference: Reference,
}

impl SequenceService {
    pub fn new(reference: Reference) -> Result<Self, TGVError> {
        Ok(Self {
            client: Client::new(),
            reference,
        })
    }

    pub async fn close(&self) -> Result<(), TGVError> {
        // Reqwest client does not need to be closed.
        Ok(())
    }

    pub async fn query_sequence(&self, region: &Region) -> Result<Sequence, TGVError> {
        let url = self
            .get_api_url(&region.contig, region.start, region.end)
            .unwrap();

        let response: UcscResponse = self.client.get(&url).send().await?.json().await?;

        Ok(Sequence {
            start: region.start,
            sequence: response.dna,
            contig: region.contig.clone(),
        })
    }

    /// start / end: 1-based, inclusive.
    fn get_api_url(&self, contig: &Contig, start: usize, end: usize) -> Result<String, TGVError> {
        match self.reference {
            Reference::Hg19 | Reference::Hg38 | Reference::UcscGenome(_) => Ok(format!(
                "https://api.genome.ucsc.edu/getData/sequence?genome={};chrom={};start={};end={}",
                self.reference.to_string(),
                contig.name,
                start - 1, // start is 0-based, inclusive.
                end
            )),
            _ => Err(TGVError::IOError("Unsupported reference".to_string())),
        }
    }
}
