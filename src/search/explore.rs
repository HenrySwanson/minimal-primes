use crate::data_structures::{Frontier, TreeTracer};
use crate::search::SearchNode;

pub trait Explore {
    fn start(node: SearchNode) -> Self;

    fn explore_one_level(&mut self, f: impl FnMut(SearchNode) -> Vec<SearchNode>) -> bool;

    fn iter(&self) -> impl Iterator<Item = &SearchNode>;

    fn min_weight(&self) -> Option<usize>;

    // some defaults

    fn len(&self) -> usize {
        self.iter().count()
    }

    fn all_simple(&self) -> bool {
        self.iter()
            .all(|node| matches!(node, SearchNode::Simple(_)))
    }

    fn print_tree_to_stdout(&self);
}

impl Explore for Frontier<SearchNode> {
    fn start(node: SearchNode) -> Self {
        Frontier::new(node)
    }

    fn explore_one_level(&mut self, f: impl FnMut(SearchNode) -> Vec<SearchNode>) -> bool {
        self.explore_one_level(f)
    }

    fn iter(&self) -> impl Iterator<Item = &SearchNode> {
        self.iter()
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn min_weight(&self) -> Option<usize> {
        self.min_weight()
    }

    fn print_tree_to_stdout(&self) {
        self.print_tree_to_stdout();
    }
}

impl Explore for TreeTracer<SearchNode> {
    fn start(node: SearchNode) -> Self {
        TreeTracer::new(node)
    }

    fn explore_one_level(&mut self, f: impl FnMut(SearchNode) -> Vec<SearchNode>) -> bool {
        self.explore_one_level(f)
    }

    fn iter(&self) -> impl Iterator<Item = &SearchNode> {
        self.iter()
    }

    fn min_weight(&self) -> Option<usize> {
        self.min_weight()
    }

    fn print_tree_to_stdout(&self) {
        self.print_tree_to_stdout();
    }
}
