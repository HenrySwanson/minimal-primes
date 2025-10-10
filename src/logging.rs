use log::{Metadata, Record};

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

/// Used to log the tree-like structure of how the search space is explored.
pub enum Tracer {
    // The tree and the current index we're looking at
    Real(AppendTree<String>, AppendTreeNodeID),
    // we need to get an id sometimes, but we can't make children or items.
    // using an empty type forces this :)
    Dummy(AppendTree<Never>),
}

/// A n-ary tree containing `T`s, where each node contains a `T` and an ordered
/// list of children, which can be more `T`s or more nodes (potentially
/// interleaved).
///
/// This allows structures like:
/// ```
/// Node(T1)   (mix of Ts and child nodes)
/// - T2
/// - T3
/// - Node(T4) (no children)
/// - Node(T5) (only Ts)
///   - T6
///   - T7
///   - T8
/// - Node(T9) (lots of T-less nodes)
///   - Node(T10)
///   - Node(T11)
/// - T12
/// ```
///
/// Used for tracking the evolution of the search tree as we look for minimal
/// primes.
#[derive(Debug, Clone)]
pub struct AppendTree<T> {
    nodes: Vec<AppendTreeNode<T>>,
}

/// Identifies nodes in a [AppendTree].
#[derive(Debug, Clone, Copy)]
pub struct AppendTreeNodeID(usize);

/// Internal nodes in an [AppendTree].
#[derive(Debug, Clone)]
struct AppendTreeNode<T> {
    contents: Vec<Content<T>>,
    // TODO: parent?
}

/// Children of an [AppendTreeNode]
#[derive(Debug, Clone)]
enum Content<T> {
    Item(T),
    Child { tag: T, idx: usize },
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

impl<T> AppendTree<T> {
    /// Creates a tree with a single root node.
    pub fn new() -> Self {
        let root = AppendTreeNode::new();
        Self { nodes: vec![root] }
    }

    /// Returns the root of this tree.
    pub fn root(&self) -> AppendTreeNodeID {
        AppendTreeNodeID(0)
    }

    /// Appends a new item to the given subtree. If the node ID doesn't exist,
    /// returns the item in the `Err` variant.
    pub fn append(&mut self, node_id: AppendTreeNodeID, item: T) -> Result<(), T> {
        match self.nodes.get_mut(node_id.0) {
            Some(node) => {
                node.contents.push(Content::Item(item));
                Ok(())
            }
            None => Err(item),
        }
    }

    /// Appends a new child node to the given subtree, with `tag` as the first
    /// item. If the node ID doesn't exist, returns the item in the `Err` variant.
    pub fn make_child(&mut self, node_id: AppendTreeNodeID, tag: T) -> Result<AppendTreeNodeID, T> {
        // gotta get this before the mutable borrow begins
        let num_nodes = self.nodes.len();

        match self.nodes.get_mut(node_id.0) {
            Some(node) => {
                // Add a new node to the whole tree, and push it into
                // this node's children.
                let child_idx = num_nodes;
                node.contents.push(Content::Child {
                    tag,
                    idx: child_idx,
                });
                self.nodes.push(AppendTreeNode::new());
                Ok(AppendTreeNodeID(child_idx))
            }
            None => Err(tag),
        }
    }
}

// is this principled? no. does it work well? yeah
impl<T: std::fmt::Display> AppendTree<T> {
    pub fn pretty_print_to_stdout(&self) {
        self.pretty_print_helper(0, 0);
    }

    fn pretty_print_helper(&self, node_idx: usize, indent: usize) {
        let node = &self.nodes[node_idx];

        for content in &node.contents {
            match content {
                Content::Item(t) => println!("{:indent$}{}", "", t, indent = indent * 2),
                Content::Child {
                    tag,
                    idx: child_idx,
                } => {
                    println!("{:indent$}{}", "", tag, indent = indent * 2);
                    self.pretty_print_helper(*child_idx, indent + 1);
                }
            }
        }
    }
}

impl<T> AppendTreeNode<T> {
    pub fn new() -> Self {
        Self { contents: vec![] }
    }
}
