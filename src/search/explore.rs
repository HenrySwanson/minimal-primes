use crate::data_structures::WeightedVec;

pub struct Frontier<T> {
    /// maps weight to nodes; used to ensure we're exploring the
    /// search space in (non-strictly) increasing order.
    by_weight: WeightedVec<T>,
}

pub trait Weight {
    fn weight(&self) -> usize;
}

impl<T: Weight> Frontier<T> {
    pub fn start(node: T) -> Self {
        let mut ret = Self {
            by_weight: WeightedVec::new(),
        };
        ret.put(node);
        ret
    }

    // TODO: return some richer type from the closure?
    pub fn explore_next(&mut self, f: impl FnOnce(T) -> Vec<T>) -> bool {
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

    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.by_weight.iter()
    }

    pub fn len(&self) -> usize {
        self.by_weight.len()
    }

    pub fn min_weight(&self) -> Option<usize> {
        self.by_weight.min_weight()
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
