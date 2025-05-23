use crate::{
    alignment::{AlignedRead, Alignment},
    error::TGVError,
    rendering::colors,
    window::{OnScreenCoordinate, ViewingWindow},
};
use ratatui::{buffer::Buffer, layout::Rect, style::Style};
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::cigar::op::Op;
use std::io;

/// Render an alignment on the alignment area.
pub fn render_alignment(
    area: &Rect,
    buf: &mut Buffer,
    window: &ViewingWindow,
    alignment: &Alignment,
) -> Result<(), TGVError> {
    // This iterates through all cached reads and re-calculates coordinates for each movement.
    // Consider improvement.
    for read in alignment.reads.iter() {
        for (x, y, onscreen_string, style) in get_read_rendering_info(read, window, area)? {
            buf.set_string(x as u16 + area.x, y as u16 + area.y, onscreen_string, style);
        }
    }
    Ok(())
}

fn get_read_rendering_info(
    read: &AlignedRead,
    viewing_window: &ViewingWindow,
    area: &Rect,
) -> io::Result<Vec<(usize, usize, String, Style)>> {
    let mut output = Vec::new();
    let cigar_segments = get_cigar_segments(read)?;
    let n_cigar_segments = cigar_segments.iter().len();

    let onscreen_y = match viewing_window.onscreen_y_coordinate(read.y, area) {
        OnScreenCoordinate::OnScreen(y_start) => y_start,
        _ => return Ok(vec![]),
    };

    for (i_cigar_segment, (start_coord, end_coord, style)) in cigar_segments.iter().enumerate() {
        if let Some((x, length)) = OnScreenCoordinate::onscreen_start_and_length(
            &viewing_window.onscreen_x_coordinate(*start_coord, area),
            &viewing_window.onscreen_x_coordinate(*end_coord, area),
            area,
        ) {
            output.push((
                x,
                onscreen_y,
                get_segment_string(length, {
                    if i_cigar_segment == 0 {
                        Some(true)
                    } else if i_cigar_segment == n_cigar_segments - 1 {
                        Some(false)
                    } else {
                        None
                    }
                }),
                *style,
            ));
        }
    }
    
    Ok(output)
}

fn get_segment_string(length: usize, is_reverse: Option<bool>) -> String {
    match is_reverse {
        Some(true) => (0..length)
            .map(|i| if i == 0 { "<" } else { "-" })
            .collect::<String>(),
        Some(false) => (0..length)
            .map(|i| if i == length - 1 { ">" } else { "-" })
            .collect::<String>(),
        None => "-".repeat(length),
    }
}

/// Render a read as sections of styled texts
/// See: https://samtools.github.io/hts-specs/SAMv1.pdf
fn get_cigar_segments(read: &AlignedRead) -> io::Result<Vec<(usize, usize, Style)>> {
    let mut reference_pivot: usize = read.start; // used in the output
    let mut query_pivot: usize = 0; // # bases relative to the softclip start.

    let mut output = Vec::new();

    for op in read.read.cigar().iter() {
        let op = op?;
        if let Kind::SoftClip = op.kind() {
            for i_base in query_pivot..query_pivot + op.len() {
                let base_coord_is_valid = reference_pivot + i_base >= 1 + read.leading_softclips;
                if base_coord_is_valid {
                    let abs_start = reference_pivot + i_base - read.leading_softclips;

                    let base = read.read.sequence().get(i_base);
                    let base_color = if let Some(base) = base {
                        match base {
                            b'A' => colors::SOFTCLIP_A,
                            b'C' => colors::SOFTCLIP_C,
                            b'G' => colors::SOFTCLIP_G,
                            b'T' => colors::SOFTCLIP_T,
                            _ => colors::SOFTCLIP_N,
                        }
                    } else {
                        colors::SOFTCLIP_N
                    };

                    output.push((abs_start, abs_start, Style::default().bg(base_color)));
                }
            }
        }

        if consumes_reference(&op) {
            output.push((
                reference_pivot,
                reference_pivot + op.len() as usize - 1_usize,
                get_cigar_style(&op),
            ));
            reference_pivot += op.len() as usize;
            // Note that softclip does not consume query and is handled above.
        }

        if consumes_query(&op) {
            query_pivot += op.len() as usize;
        }
    }

    Ok(output)
}

/// Whether the cigar operation consumes reference.
/// Yes: M/D/N/=/X
/// No: I/S/H/P
/// See: https://samtools.github.io/hts-specs/SAMv1.pdf
fn consumes_reference(op: &Op) -> bool {
    match op.kind() {
        Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch => true,
        Kind::SoftClip | Kind::Insertion | Kind::HardClip | Kind::Pad => false
    }
}

/// Whether the cigar operation consumes query.
/// Yes: M/I/S/=/X
/// No: D/N/H/P
fn consumes_query(op: &Op) -> bool {
    match op.kind() {
        Kind::Match | Kind::Insertion | Kind::SoftClip | Kind::SequenceMatch | Kind::SequenceMismatch => true,
        Kind::Deletion | Kind::Skip | Kind::HardClip | Kind::Pad => false,
    }
}

/// Only labels that consumes reference are display onscreen.
fn get_cigar_style(op: &Op) -> Style {
    match op.kind() {
        Kind::Match | Kind::SequenceMatch => Style::default().bg(colors::MATCH_COLOR),
        // By SAM spec, M can also be mismatch. TODO: think about this in the future.
        Kind::SequenceMismatch => Style::default().bg(colors::MISMATCH_COLOR),
        _ => Style::default(),
    }
}
