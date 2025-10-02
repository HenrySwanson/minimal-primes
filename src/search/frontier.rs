use std::ops::ControlFlow;

use crate::data_structures::WeightedVec;

/// Holds the families we haven't explored yet, arranged by increasing
/// weight. Helps guarantee we've discovered all minimal primes of length N
/// before examining strings of length N+1.
///
/// Holds [super::SearchNode]s for our purposes, but in theory can hold
/// anything implementing [Weight].
pub struct Frontier<T> {
    /// maps weight to nodes; used to ensure we're exploring the
    /// search space in (non-strictly) increasing order.
    by_weight: WeightedVec<T>,
}

/// Helper trait to make [Frontier] work.
pub trait Weight {
    fn weight(&self) -> usize;
}

impl<T: Weight> Frontier<T> {
    /// Creates a new frontier with exactly one element.
    pub fn start(node: T) -> Self {
        let mut ret = Self {
            by_weight: WeightedVec::new(),
        };
        ret.put(node);
        ret
    }

    /// Pops an element of least weight from the frontier, passes it to the closure,
    /// and inserts the output into the frontier.
    ///
    /// Returns [ControlFlow::Break] if there is nothing left in the frontier.
    pub fn explore_next(&mut self, f: impl FnOnce(T) -> Vec<T>) -> ControlFlow<()> {
        // TODO: return some richer type from the closure?

        // Pop out an element of least weight
        let layer = match self.by_weight.find_first_non_empty_layer_mut() {
            Some(layer) => layer,
            None => return ControlFlow::Break(()),
        };

        let node = layer.pop_front().expect("non-empty layer");
        for child in f(node) {
            self.put(child);
        }

        ControlFlow::Continue(())
    }

    /// Iterates through every item in the frontier.
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.by_weight.iter()
    }

    /// Returns the number of elements in the frontier.
    pub fn len(&self) -> usize {
        self.by_weight.len()
    }

    /// Returns the minimum weight across all items in the frontier,
    /// or None if the frontier is empty.
    pub fn min_weight(&self) -> Option<usize> {
        self.by_weight.min_weight()
    }
    
    /// Inserts a new item into the frontier.
    fn put(&mut self, node: T) {
        let weight = node.weight();
        self.by_weight.put(node, weight);
    }
}
