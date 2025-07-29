use std::fmt::Display;

use crate::data_structures::WeightedVec;
use crate::search::{NodeType, SearchNode};

pub trait Explore {
    fn start(node: SearchNode) -> Self;

    /// Pops the next node, and applies the given function to it,
    /// inserting the returned nodes into itself.
    ///
    /// Returned nodes must have weight at least as large as the
    /// popped node.
    ///
    /// Returns false if this structure is empty and there are no
    /// nodes to explore.
    fn explore_next(&mut self, f: impl FnOnce(SearchNode) -> Vec<SearchNode>) -> bool;

    fn iter(&self) -> impl Iterator<Item = &SearchNode>;

    fn min_weight(&self) -> Option<usize>;

    // some defaults

    fn len(&self) -> usize {
        self.iter().count()
    }

    fn all_simple(&self) -> bool {
        self.iter()
            .all(|node| matches!(node.family, NodeType::Simple(_)))
    }

    fn print_tree_to_stdout(&self);
}

pub struct Frontier<T> {
    /// maps weight to nodes; used to ensure we're exploring the
    /// search space in (non-strictly) increasing order.
    by_weight: WeightedVec<T>,
}

pub struct TreeTracer<T> {
    nodes: Vec<TreeTracerNode<T>>,
    unexplored: WeightedVec<usize>,
}

struct TreeTracerNode<T> {
    value: T,
    reason: Option<String>,
    children: Vec<usize>,
}

pub trait Weight {
    fn weight(&self) -> usize;
}

impl Explore for Frontier<SearchNode> {
    fn start(node: SearchNode) -> Self {
        let mut ret = Self {
            by_weight: WeightedVec::new(),
        };
        ret.put(node);
        ret
    }

    // TODO: return some richer type from the closure?
    fn explore_next(&mut self, f: impl FnOnce(SearchNode) -> Vec<SearchNode>) -> bool {
        // Pop out an element of least weight
        let layer = match self.by_weight.find_first_non_empty_layer_mut() {
            Some(layer) => layer,
            None => return false,
        };

        let node = layer.pop_front().expect("non-empty layer");
        for child in f(node) {
            self.put(child);
        }
        true
    }

    fn iter(&self) -> impl Iterator<Item = &SearchNode> {
        self.by_weight.iter()
    }

    fn len(&self) -> usize {
        self.by_weight.len()
    }

    fn min_weight(&self) -> Option<usize> {
        self.by_weight.min_weight()
    }

    fn print_tree_to_stdout(&self) {
        // do nothing
    }
}

impl<T> Frontier<T>
where
    T: Weight,
{
    fn put(&mut self, node: T) {
        let weight = node.weight();
        self.by_weight.put(node, weight);
    }
}

impl Explore for TreeTracer<SearchNode> {
    fn start(node: SearchNode) -> Self {
        let mut ret = Self {
            nodes: vec![],
            unexplored: WeightedVec::new(),
        };
        ret.put(node);
        ret
    }

    fn explore_next(&mut self, f: impl FnOnce(SearchNode) -> Vec<SearchNode>) -> bool {
        // Pop out an element of least weight
        let layer = match self.unexplored.find_first_non_empty_layer_mut() {
            Some(layer) => layer,
            None => return false,
        };

        let node_idx = layer.pop_front().expect("non-empty layer");
        let node = &self.nodes[node_idx];
        let mut child_idxs = vec![];
        for child in f(node.value.clone()) {
            let child_idx = self.nodes.len();
            child_idxs.push(child_idx);

            self.put(child);
        }

        self.nodes[node_idx].children = child_idxs;
        self.nodes[node_idx].reason = Some(String::new());
        true
    }

    fn iter(&self) -> impl Iterator<Item = &SearchNode> {
        self.unexplored.iter().map(|idx| &self.nodes[*idx].value)
    }

    fn len(&self) -> usize {
        self.unexplored.len()
    }

    fn min_weight(&self) -> Option<usize> {
        self.unexplored.min_weight()
    }

    fn print_tree_to_stdout(&self) {
        println!("---- SEARCH TREE ----");
        self.print_tree_to_stdout_helper(0, 0);
    }
}

impl<T> TreeTracer<T> {
    fn put(&mut self, node: T)
    where
        T: Weight,
    {
        let weight = node.weight();
        let idx = self.nodes.len();
        self.nodes.push(TreeTracerNode {
            value: node,
            reason: None,
            children: vec![],
        });
        self.unexplored.put(idx, weight);
    }

    fn print_tree_to_stdout_helper(&self, node_idx: usize, indent: usize)
    where
        T: Display,
    {
        let node = &self.nodes[node_idx];

        match &node.reason {
            Some(reason) => {
                println!(
                    "{:indent$}{}  {}",
                    "",
                    node.value,
                    reason,
                    indent = indent * 2,
                );
            }
            None => {
                println!("{:indent$}{}", "", node.value, indent = indent * 2);
            }
        }
        for child_idx in &node.children {
            self.print_tree_to_stdout_helper(*child_idx, indent + 1);
        }
    }
}
