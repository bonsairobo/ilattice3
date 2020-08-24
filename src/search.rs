// This code is copied from the pathfinding crate and modified to support a specific use case.
// Licensed under dual MIT / Apache 2.0 at the time of copying.

use indexmap::map::Entry::{Occupied, Vacant};
use indexmap::IndexMap;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::hash::Hash;

pub fn greedy_best_first<N, C, FN, IN, FH, FS>(
    start: &N,
    mut successors: FN,
    mut heuristic: FH,
    mut success: FS,
    max_iterations: usize,
) -> (bool, Vec<N>)
where
    N: Eq + Hash + Clone,
    C: Ord + Copy,
    FN: FnMut(&N) -> IN,
    IN: IntoIterator<Item = N>,
    FH: FnMut(&N) -> C,
    FS: FnMut(&N) -> bool,
{
    let h_start = heuristic(start);
    let mut best_heuristic_so_far = h_start;
    let mut best_index_so_far = 0;

    let mut to_see = BinaryHeap::new();
    to_see.push(HeuristicCostHolder {
        estimated_cost: h_start,
        index: 0,
    });
    let mut parents: IndexMap<N, usize> = IndexMap::new();
    parents.insert(start.clone(), usize::max_value());
    let mut num_iters = 0;
    while let Some(HeuristicCostHolder { index, .. }) = to_see.pop() {
        let successors = {
            let (node, _) = parents.get_index(index).unwrap();
            if success(node) {
                let path = reverse_path(&parents, index);
                return (true, path);
            }

            successors(node)
        };

        for successor in successors {
            let h; // heuristic(&successor)
            let n; // index for successor
            match parents.entry(successor) {
                Vacant(e) => {
                    h = heuristic(e.key());
                    n = e.index();
                    e.insert(index);
                }
                Occupied(_) => {
                    continue;
                }
            }

            to_see.push(HeuristicCostHolder {
                estimated_cost: h,
                index: n,
            });

            if h < best_heuristic_so_far {
                best_heuristic_so_far = h;
                best_index_so_far = n;
            }
        }

        num_iters += 1;
        if num_iters == max_iterations {
            break;
        }
    }

    let path = reverse_path(&parents, best_index_so_far);

    (false, path)
}

fn reverse_path<N>(parents: &IndexMap<N, usize>, start: usize) -> Vec<N>
where
    N: Eq + Hash + Clone,
{
    let path = itertools::unfold(start, |i| {
        parents.get_index(*i).map(|(node, index)| {
            *i = *index;

            node
        })
    })
    .collect::<Vec<&N>>();

    path.into_iter().rev().cloned().collect()
}

struct HeuristicCostHolder<K> {
    estimated_cost: K,
    index: usize,
}

impl<K: PartialEq> PartialEq for HeuristicCostHolder<K> {
    fn eq(&self, other: &Self) -> bool {
        self.estimated_cost.eq(&other.estimated_cost)
    }
}

impl<K: PartialEq> Eq for HeuristicCostHolder<K> {}

impl<K: Ord> PartialOrd for HeuristicCostHolder<K> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<K: Ord> Ord for HeuristicCostHolder<K> {
    fn cmp(&self, other: &Self) -> Ordering {
        other.estimated_cost.cmp(&self.estimated_cost)
    }
}
