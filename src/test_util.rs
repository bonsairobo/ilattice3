use std::collections::HashSet;
use std::fmt::Debug;
use std::hash::Hash;
use std::iter::FromIterator;

pub fn assert_elements_eq<T: Clone + Debug + Eq + Hash>(v1: &Vec<T>, v2: &Vec<T>) {
    let set1: HashSet<T> = HashSet::from_iter(v1.iter().cloned());
    let set2: HashSet<T> = HashSet::from_iter(v2.iter().cloned());
    assert_eq!(set1, set2);
}
