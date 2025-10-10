use std::collections::VecDeque;
use std::ops::ControlFlow;

/// Holds the families we haven't explored yet, arranged by increasing
/// weight. Helps guarantee we've discovered all minimal primes of length N
/// before examining strings of length N+1.
///
/// Holds [super::SearchNode]s for our purposes, but in theory can hold
/// anything implementing [Weight].
pub struct Frontier<T> {
    /// maps weight to nodes; an element with weight i is in the ith deque.
    /// used to ensure we're exploring the search space in (non-strictly)
    /// increasing order
    elements: Vec<VecDeque<T>>,
    /// maps weight to nodes, the 'ratchet' that enforces that we can't backtrack
    /// to an element of lower weight.
    min_allowed_weight: usize,
}

/// Helper trait to make [Frontier] work.
pub trait Weight {
    fn weight(&self) -> usize;
}

impl<T: Weight> Frontier<T> {
    /// Creates a new frontier with exactly one element.
    pub fn start(node: T) -> Self {
        let mut ret = Self {
            elements: vec![],
            min_allowed_weight: 0,
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
        let layer = match self.find_first_non_empty_layer_mut() {
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
        self.elements.iter().flatten()
    }

    /// Returns the number of elements in the frontier.
    pub fn len(&self) -> usize {
        self.elements.iter().map(|layer| layer.len()).sum()
    }

    /// Returns the minimum weight across all items in the frontier,
    /// or None if the frontier is empty.
    pub fn min_weight(&self) -> Option<usize> {
        self.elements.iter().position(|layer| !layer.is_empty())
    }

    /// Inserts a new item into the frontier.
    fn put(&mut self, item: T) {
        let weight = item.weight();
        debug_assert!(weight >= self.min_allowed_weight);

        loop {
            // Keep appending empty deques until we reach the right weight
            match self.elements.get_mut(weight) {
                Some(layer) => {
                    layer.push_back(item);
                    break;
                }
                None => self.elements.push(VecDeque::new()),
            }
        }
    }

    fn find_first_non_empty_layer_mut(&mut self) -> Option<&mut VecDeque<T>> {
        for (i, layer) in self.elements.iter_mut().enumerate() {
            if i < self.min_allowed_weight {
                assert!(
                    layer.is_empty(),
                    "all layers below min_allowed_weight must be empty"
                );
            }

            if !layer.is_empty() {
                return Some(layer);
            }
        }

        None
    }
}
