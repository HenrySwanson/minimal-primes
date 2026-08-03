//! Reconstructs and pretty-prints the search tree from a `--trace` JSONL
//! file, in natural tree order rather than the (weight-ordered) order it was
//! actually explored in.
//!
//! Because a node is recorded in the trace *after* it's been explored, if the
//! run is cut off partway through, any children that were created but not yet
//! explored, they will not show up in the tree. That's fine enough, they can
//! be inferred from the parent's event if necessary.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

use itertools::Itertools;

use crate::search::{CompositeReason, ExploreEvent, SplitDirection, TraceRecord};

/// The root of every search tree is always the first node ever allocated.
const ROOT_ID: u64 = 0;

pub fn print_tree(path: &Path) -> io::Result<()> {
    let records = read_trace(path)?;

    if !records.0.contains_key(&ROOT_ID) {
        println!("No record for the root node (id {ROOT_ID}) -- is this a complete trace file?");
        return Ok(());
    }

    // Grab the stdout lock up front -- we're gonna be doing a lot of printing.
    let stdout = io::stdout();
    print_subtree(ROOT_ID, &records, &mut stdout.lock())
}

/// Maps a node ID to its record and its children.
struct TraceRecords(HashMap<u64, (TraceRecord, Vec<u64>)>);

/// Reads every record out of a trace file, and indexes them by id, as well
/// as by their parent's id (so we can find a node's children).
fn read_trace(path: &Path) -> io::Result<TraceRecords> {
    let file = File::open(path)?;

    let mut records: HashMap<u64, (TraceRecord, Vec<u64>)> = HashMap::new();

    for line in BufReader::new(file).lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }

        let record: TraceRecord = serde_json::from_str(&line)
            .unwrap_or_else(|e| panic!("failed to parse trace record: {e}\n  {line}"));

        // If we've got a parent (almost certainly), put our node ID into the parents'
        // map entry.
        if let Some(parent_id) = record.parent_id {
            let Some(parent_value) = records.get_mut(&parent_id) else {
                panic!(
                    "Record {} has parent {}, which was not yet seen",
                    record.node_id, parent_id
                );
            };
            parent_value.1.push(record.node_id);
        }

        // Then put ourselves in the map
        records.insert(record.node_id, (record, vec![]));
    }

    Ok(TraceRecords(records))
}

/// Walks the tree rooted at `id` in natural (pre-order) order, printing one
/// line per node with unicode box-drawing connectors.
///
/// This is iterative, rather than the more natural recursion, because some
/// branches of the search tree (e.g. a slowly-incrementing simple family)
/// can be hundreds of thousands of nodes deep, which would blow the native
/// call stack.
fn print_subtree(root: u64, records: &TraceRecords, out: &mut impl Write) -> io::Result<()> {
    // We'll maintain a stack of triples:
    // - the node id
    // - the prefix to print before its own line
    // - the prefix to print before its children's lines
    //
    // Those last two can be different! When printing a child, you need to know
    // whether it's the last of its siblings (to get the right-angle connector instead
    // of the T-connector), but also whether its parent is the last of its siblings
    // and same with its grandparents, and so on. So we just track them separately.
    let mut stack = vec![(root, String::new(), String::new())];

    while let Some((id, line_prefix, child_prefix)) = stack.pop() {
        let (record, children) = &records.0[&id];

        writeln!(out, "{line_prefix}{}", describe_node(record))?;

        // Push in reverse, so the first child ends up on top of the
        // stack and gets printed first.
        for (i, child_id) in children.iter().rev().copied().enumerate() {
            let is_last = i == 0;
            let connector = if is_last { "└── " } else { "├── " };
            let continuation = if is_last { "    " } else { "│   " };
            stack.push((
                child_id,
                format!("{child_prefix}{connector}"),
                format!("{child_prefix}{continuation}"),
            ));
        }
    }

    Ok(())
}

fn describe_node(record: &TraceRecord) -> String {
    let event = describe_event(&record.event);
    match &record.reduced {
        Some(reduced) => format!("{}: reduced to {reduced}, {event}", record.family),
        None => format!("{}: {event}", record.family),
    }
}

fn describe_event(event: &ExploreEvent) -> String {
    match event {
        ExploreEvent::ContainsPrime(p) => format!("contains known prime {p}"),
        ExploreEvent::IsNewPrime => "is a new prime!".to_string(),
        ExploreEvent::NoCoresRemaining => "reduced to a trivial string".to_string(),
        ExploreEvent::DetectedComposite(reason) => {
            format!("detected composite ({})", describe_composite_reason(reason))
        }
        ExploreEvent::Simplified => "simplified to a single-core family".to_string(),
        ExploreEvent::SplitOnLimitedDigit { core_idx, digit, n } => {
            format!("split on limited digit ({digit} can't repeat {n} times in core {core_idx})")
        }
        ExploreEvent::SplitOnIncompatibleDifferentCores {
            core_i,
            core_j,
            a,
            b,
        } => format!(
            "split on incompatible digits ({a} in core {core_i}, {b} in core {core_j} can't co-occur)"
        ),
        ExploreEvent::SplitOnIncompatibleSameCore {
            core_idx,
            first,
            second,
            reverse_also_forbidden,
        } => {
            if *reverse_also_forbidden {
                format!(
                    "split on incompatible digits ({first} and {second} can't co-occur in core {core_idx})"
                )
            } else {
                format!(
                    "split on incompatible digits ({first}{second} can't appear in core {core_idx})"
                )
            }
        }
        ExploreEvent::SplitOnForbiddenSandwich { core_idx, a, b } => {
            format!("split on forbidden sandwich ({a}{b}{a} in core {core_idx})")
        }
        ExploreEvent::SplitOnNecessaryDigit { core_idx, digit } => {
            format!("split on necessary digit ({digit} required in core {core_idx})")
        }
        ExploreEvent::SplitGenerically {
            core_idx,
            direction,
        } => {
            let direction = match direction {
                SplitDirection::Left => "left",
                SplitDirection::Right => "right",
            };
            format!("split generically {direction} on core {core_idx}")
        }
        ExploreEvent::IncrementedRepeat => "incremented repeat count".to_string(),
    }
}

fn describe_composite_reason(reason: &CompositeReason) -> String {
    match reason {
        CompositeReason::SharesFactorWithBase(p) => format!("shares factor {p} with base"),
        CompositeReason::CommonFactor(p) => format!("common factor of {p}"),
        CompositeReason::PeriodicFactors(factors) => {
            format!("periodic factors {}", factors.iter().format(", "))
        }
        CompositeReason::LocalAlternatingFactors {
            core_idx,
            even_factor,
            odd_factor,
        } => format!("alternating factors {even_factor}/{odd_factor} in core {core_idx}"),
        CompositeReason::GlobalAlternatingFactors {
            even_factor,
            odd_factor,
        } => format!("alternating factors {even_factor}/{odd_factor}"),
        CompositeReason::NeverCoprimeTo30 => "never coprime to 30".to_string(),
        CompositeReason::FactorsAlgebraically => "factors algebraically".to_string(),
    }
}
