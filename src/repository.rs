use crate::{
    alignment::{Alignment, AlignmentBuilder},
    contig::Contig,
    error::TGVError,
    helpers::is_url,
    reference::Reference,
    region::Region,
    sequence::Sequence,
    settings::{BackendType, Settings},
    track_service::{TrackService, TrackServiceEnum, UcscDbTrackService},
};
use noodles::bam::io::indexed_reader;
use noodles::bam::io::reader;

use reqwest::Client;
use noodles::bam::{self, io::IndexedReader};
use serde::Deserialize;
use std::path::Path;
use url::Url;
use noodles::bam::bai;
use noodles::sam::Header;
use noodles::core::region;
use noodles::core::Position;

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
enum RemoteSource {
    S3,
    HTTP,
    GS,
}

impl RemoteSource {
    fn from(path: &String) -> Result<Self, TGVError> {
        if path.starts_with("s3://") {
            Ok(Self::S3)
        } else if path.starts_with("http://") || path.starts_with("https://") {
            Ok(Self::HTTP)
        } else if path.starts_with("gss://") {
            Ok(Self::GS)
        } else {
            Err(TGVError::ValueError(format!(
                "Unsupported remote path {}. Only S3, HTTP/HTTPS, and GS are supported.",
                path
            )))
        }
    }
}

pub trait AlignmentRepository {
    fn read_alignment(&self, region: &Region) -> Result<Alignment, TGVError>;

    fn read_header(&self) -> Result<Vec<(String, Option<usize>)>, TGVError>;
}

#[derive(Debug)]
pub struct BamRepository {
    bam_path: String,
    bai_path: Option<String>,
}

impl BamRepository {
    fn new(bam_path: String, bai_path: Option<String>) -> Result<Self, TGVError> {
        if is_url(&bam_path) {
            return Err(TGVError::IOError(format!(
                "{} is a remote path. Use RemoteBamRepository for remote BAM IO",
                bam_path
            )));
        }

        if !Path::new(&bam_path).exists() {
            return Err(TGVError::IOError(format!(
                "BAM file {} not found",
                bam_path
            )));
        }

        match &bai_path {
            Some(bai_path) => {
                if !Path::new(bai_path).exists() {
                    return Err(TGVError::IOError(format!(
                        "BAM index file {} not found. Only indexed BAM files are supported.",
                        bai_path
                    )));
                }
            }
            None => {
                if !Path::new(&format!("{}.bai", bam_path)).exists() {
                    return Err(TGVError::IOError(format!(
                        "BAM index file {}.bai not found. Only indexed BAM files are supported.",
                        bam_path
                    )));
                }
            }
        }

        Ok(Self { bam_path, bai_path })
    }
}

impl AlignmentRepository for BamRepository {
    fn read_alignment(&self, region: &Region) -> Result<Alignment, TGVError> {
        let mut reader = match self.bai_path.as_ref() {
            Some(bai_path) => {
                let index = bai::fs::read(bai_path)?;
                indexed_reader::Builder::default().set_index(index).build_from_path(self.bam_path.clone())?
            }
            None => {
                indexed_reader::Builder::default().build_from_path(self.bam_path.clone())?
            }
        };

        let header = reader.read_header()?;
        let noodles_region = region::Region::new(
            region.contig.name.to_string(),
            Position::new(region.start).ok_or_else(|| TGVError::ValueError("invalid position".to_string()))?..=Position::new(region.end).ok_or_else(|| TGVError::ValueError("invalid position".to_string()))?,
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
        let mut reader = match self.bai_path.as_ref() {
            Some(bai_path) => {
                let index = bai::fs::read(bai_path)?;
                indexed_reader::Builder::default().set_index(index).build_from_path(self.bam_path.clone())?
            }
            None => {
                indexed_reader::Builder::default().build_from_path(self.bam_path.clone())?
            }
        };

        let header = reader.read_header()?;
        get_contig_names_and_lengths_from_header(&header)
    }
}

#[derive(Debug)]
pub struct RemoteBamRepository {
    bam_path: String,
    source: RemoteSource,
}

impl RemoteBamRepository {
    pub fn new(bam_path: &String) -> Result<Self, TGVError> {
        Ok(Self {
            bam_path: bam_path.clone(),
            source: RemoteSource::from(bam_path)?,
        })
    }
}

impl AlignmentRepository for RemoteBamRepository {
    fn read_alignment(&self, region: &Region) -> Result<Alignment, TGVError> {
        let mut bam = IndexedReader::from_url(
            &Url::parse(&self.bam_path).map_err(|e| TGVError::IOError(e.to_string()))?,
        )?;

        let header = bam::Header::from_template(bam.header());

        let query_contig_string = get_query_contig_string(&header, region)?;
        bam.fetch((
            &query_contig_string,
            region.start as i32 - 1,
            region.end as i32,
        ))
        .map_err(|e| TGVError::IOError(e.to_string()))?;

        let mut alignment_builder = AlignmentBuilder::new()?;

        for record in bam.records() {
            let read = record.map_err(|e| TGVError::IOError(e.to_string()))?;
            alignment_builder.add_read(read)?;
        }

        alignment_builder.region(region)?.build()
    }

    fn read_header(&self) -> Result<Vec<(String, Option<usize>)>, TGVError> {
        let bam = IndexedReader::from_url(
            &Url::parse(&self.bam_path).map_err(|e| TGVError::IOError(e.to_string()))?,
        )?;

        let header = bam::Header::from_template(bam.header());
        get_contig_names_and_lengths_from_header(&header)
    }
}

// fn is_remote_path {
//     IndexedReader::from_url(
//         &Url::parse(bam_path).map_err(|e| TGVError::IOError(e.to_string()))?,
//     )
//     .unwrap();

// struct CRAMRepository {
//     cram_path: String,
// }

fn get_contig_names_and_lengths_from_header(
    header: &Header,
) -> Result<Vec<(String, Option<usize>)>, TGVError> {
    let mut output = Vec::new();

    // header.reference_sequences()

    for (_key, reference_sequence) in header.reference_sequences() {
        for record in reference_sequence.other_fields() {
            // match record.0 {
            //     Some(Standard
            // }

            if record.contains_key("SN") {
                let contig_name = record["SN"].to_string();
                let contig_length = if record.contains_key("LN") {
                    record["LN"].to_string().parse::<usize>().ok()
                } else {
                    None
                };

                output.push((contig_name, contig_length))
            }
        }
    }

    Ok(output)
}

#[derive(Debug)]
pub enum AlignmentRepositoryEnum {
    None,
    Bam(BamRepository),
    RemoteBam(RemoteBamRepository),
}

impl AlignmentRepositoryEnum {
    pub fn from(settings: &Settings) -> Result<Self, TGVError> {
        if settings.bam_path.is_none() {
            return Ok(AlignmentRepositoryEnum::None);
        }

        let bam_path = settings.bam_path.clone().unwrap();

        if is_url(&bam_path) {
            return Ok(AlignmentRepositoryEnum::RemoteBam(
                RemoteBamRepository::new(&bam_path)?,
            ));
        }

        Ok(AlignmentRepositoryEnum::Bam(BamRepository::new(
            bam_path,
            settings.bai_path.clone(),
        )?))
    }

    pub fn has_alignment(&self) -> bool {
        match self {
            AlignmentRepositoryEnum::Bam(_) => true,
            AlignmentRepositoryEnum::RemoteBam(_) => true,
            AlignmentRepositoryEnum::None => false,
        }
    }
}

impl AlignmentRepository for AlignmentRepositoryEnum {
    fn read_alignment(&self, region: &Region) -> Result<Alignment, TGVError> {
        match self {
            AlignmentRepositoryEnum::Bam(repository) => repository.read_alignment(region),
            AlignmentRepositoryEnum::RemoteBam(repository) => repository.read_alignment(region),
            AlignmentRepositoryEnum::None => Err(TGVError::IOError("No alignment".to_string())),
        }
    }

    fn read_header(&self) -> Result<Vec<(String, Option<usize>)>, TGVError> {
        match self {
            AlignmentRepositoryEnum::Bam(repository) => repository.read_header(),
            AlignmentRepositoryEnum::RemoteBam(repository) => repository.read_header(),
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
