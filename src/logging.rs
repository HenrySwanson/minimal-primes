use log::{Metadata, Record};

use crate::data_structures::{AppendTree, AppendTreeNodeID};

pub struct SimpleLogger;

impl log::Log for SimpleLogger {
    fn enabled(&self, _metadata: &Metadata) -> bool {
        true
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            println!("{} - {}", record.level(), record.args());
        }
    }

    fn flush(&self) {}
}

pub enum Never {}

pub enum Tracer {
    // The tree and the current index we're looking at
    Real(AppendTree<String>, AppendTreeNodeID),
    // we need to get an id sometimes, but we can't make children or items.
    // using an empty type forces this :)
    Dummy(AppendTree<Never>),
}

impl Tracer {
    pub fn new() -> Self {
        let tree = AppendTree::new();
        let root = tree.root();
        Self::Real(tree, root)
    }

    pub fn dummy() -> Self {
        Self::Dummy(AppendTree::new())
    }

    pub fn root(&self) -> AppendTreeNodeID {
        match self {
            Tracer::Real(t, _) => t.root(),
            Tracer::Dummy(t) => t.root(),
        }
    }

    pub fn make_child(
        &mut self,
        node_id: AppendTreeNodeID,
        tag: String,
    ) -> Result<AppendTreeNodeID, String> {
        match self {
            Tracer::Real(t, _) => t.make_child(node_id, tag),
            // we're never going to log anything, so just keep returning the root
            // to satisfy the type system
            Tracer::Dummy(t) => Ok(t.root()),
        }
    }

    pub fn set_id(&mut self, node_id: AppendTreeNodeID) {
        match self {
            Tracer::Real(_, id) => *id = node_id,
            Tracer::Dummy(_) => {}
        }
    }

    pub fn log(&mut self, item: String) {
        match self {
            Tracer::Real(t, id) => {
                t.append(*id, item).expect("logging to nonexistent id");
            }
            Tracer::Dummy(_) => {}
        }
    }
}
