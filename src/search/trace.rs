use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::time::{SystemTime, UNIX_EPOCH};

use serde::Serialize;

use crate::search::context::ExploreEvent;

/// Appends one JSON-lines record per explored search-tree node to a file
/// under `results/`, when tracing is enabled. Writing is buffered; if a
/// write ever fails, we log a warning and drop that record rather than take
/// down the whole search over a broken trace file.
pub struct TraceWriter {
    writer: BufWriter<File>,
}

/// One line of the trace: what happened when we explored `node_id` (a child
/// of `parent_id`, or the root if `parent_id` is `None`), and which new
/// nodes it produced.
#[derive(Serialize)]
pub struct TraceRecord {
    pub node_id: u64,
    pub parent_id: Option<u64>,
    pub family: String,
    pub event: ExploreEvent,
}

impl TraceWriter {
    /// Creates `results/` if it doesn't already exist, and opens a fresh
    /// trace file inside it, named after the base and the current time so
    /// repeated runs don't clobber each other's traces.
    pub fn create(base: u8) -> io::Result<Self> {
        std::fs::create_dir_all("results")?;

        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock should be after 1970")
            .as_secs();
        let path = format!("results/trace-base{base}-{timestamp}.jsonl");

        let file = File::create(&path)?;
        println!("Writing trace events to {path}");
        Ok(Self {
            writer: BufWriter::new(file),
        })
    }

    pub fn record(&mut self, record: &TraceRecord) {
        if let Err(e) = serde_json::to_writer(&mut self.writer, record) {
            log::warn!("failed to write trace record: {e}");
            return;
        }
        if let Err(e) = self.writer.write_all(b"\n") {
            log::warn!("failed to write trace record newline: {e}");
        }
    }
}
