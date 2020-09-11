use core::hash::{BuildHasher, Hash};
use std::collections::{
    hash_map::{self, Entry},
    HashMap,
};

/// A cache that tracks the Least Recently Used element for next eviction. Here, "used" means read
/// or written.
///
/// Eviction does not happen inline; the user must explicitly call `evict_lru` to remove the LRU
/// element. Thus, the cache may grow unbounded unless evictions or explicit removals occur.
#[derive(Clone, Debug)]
pub struct LruCache<K, V, H> {
    store: HashMap<K, usize, H>,
    order: LruList<(K, V)>,
}

impl<K: Hash + Eq, V, H> LruCache<K, V, H> {
    pub fn with_hasher(hasher_builder: H) -> LruCache<K, V, H> {
        LruCache {
            store: HashMap::with_hasher(hasher_builder),
            order: LruList::<(K, V)>::new(),
        }
    }
}

impl<K, V, H> LruCache<K, V, H>
where
    K: Hash + Eq + Clone,
    H: BuildHasher,
{
    pub fn get_mut(&mut self, key: &K) -> std::option::Option<&mut V> {
        if let Some(&index) = self.store.get(key) {
            self.order.move_to_front(index);

            Some(&mut self.order.get_mut(index).1)
        } else {
            None
        }
    }

    pub fn get(&mut self, key: &K) -> Option<&V> {
        self.get_mut(key).map(|v| &*v)
    }

    /// Allows us to get a const reference without having `&mut self`. WARNING: This will not update
    /// the LRU order!
    pub fn get_const(&self, key: &K) -> Option<&V> {
        if let Some(&index) = self.store.get(key) {
            Some(&self.order.get(index).1)
        } else {
            None
        }
    }

    pub fn insert(&mut self, key: K, val: V) -> Option<V> {
        let Self { store, order, .. } = self;
        let mut old_val = None;
        let entry = store.entry(key.clone());
        match entry {
            hash_map::Entry::Occupied(entry) => {
                let index = *entry.get();
                order.move_to_front(index);
                old_val = order.set(index, (key, val)).map(|(_, v)| v);
            }
            hash_map::Entry::Vacant(entry) => {
                entry.insert(order.push_front(Some((key, val))));
            }
        }

        old_val
    }

    pub fn get_or_insert_with(&mut self, key: K, f: impl FnOnce() -> V) -> &mut V {
        let val = self.store.entry(key);
        let Self { order, .. } = self;

        match val {
            Entry::Occupied(occupied) => {
                let index = *occupied.get();
                order.move_to_front(index);

                &mut order.get_mut(index).1
            }
            Entry::Vacant(vacant) => {
                let key = vacant.key().clone();
                let index = *vacant.insert(order.push_front(None));
                order.set(index, (key, f()));

                &mut order.get_mut(index).1
            }
        }
    }

    pub fn remove(&mut self, k: &K) -> Option<V> {
        if let Some(index) = self.store.remove(k) {
            let (_key, value) = self.order.remove(index);

            Some(value)
        } else {
            None
        }
    }

    pub fn evict_lru(&mut self) -> Option<(K, V)> {
        if self.is_empty() {
            return None;
        }

        let (key, value) = self.order.pop_back();
        self.store.remove(&key).unwrap();

        Some((key, value))
    }

    pub fn clear(&mut self) {
        self.store.clear();
        self.order.clear();
    }

    pub fn len(&self) -> usize {
        self.store.len()
    }

    pub fn is_empty(&self) -> bool {
        self.store.is_empty()
    }

    pub fn keys(&self) -> impl Iterator<Item = &K> {
        self.store.keys()
    }

    pub fn values(&self) -> impl Iterator<Item = &V> {
        let Self { store, order, .. } = self;

        store.iter().map(move |(_k, index)| &order.get(*index).1)
    }
}

/// Doubly-linked list using Vec as storage.
#[derive(Clone, Debug)]
struct LruList<T> {
    entries: Vec<ListEntry<T>>,
}

#[derive(Clone, Debug)]
struct ListEntry<T> {
    value: Option<T>,
    next: usize,
    prev: usize,
}

/// Free and occupied cells are each linked into a cyclic list with one auxiliary cell.
/// Cell #0 is on the list of free cells, element #1 is on the list of occupied cells.
impl<T> LruList<T> {
    const FREE: usize = 0;
    const OCCUPIED: usize = 1;

    fn new() -> LruList<T> {
        let mut entries = Vec::with_capacity(2);
        entries.push(ListEntry::<T> {
            value: None,
            next: 0,
            prev: 0,
        });
        entries.push(ListEntry::<T> {
            value: None,
            next: 1,
            prev: 1,
        });

        LruList { entries }
    }

    fn unlink(&mut self, index: usize) {
        let prev = self.entries[index].prev;
        let next = self.entries[index].next;
        self.entries[prev].next = next;
        self.entries[next].prev = prev;
    }

    fn link_after(&mut self, index: usize, prev: usize) {
        let next = self.entries[prev].next;
        self.entries[index].prev = prev;
        self.entries[index].next = next;
        self.entries[prev].next = index;
        self.entries[next].prev = index;
    }

    fn move_to_front(&mut self, index: usize) {
        self.unlink(index);
        self.link_after(index, Self::OCCUPIED);
    }

    fn push_front(&mut self, value: Option<T>) -> usize {
        if self.entries[Self::FREE].next == Self::FREE {
            self.entries.push(ListEntry::<T> {
                value: None,
                next: Self::FREE,
                prev: Self::FREE,
            });
            self.entries[Self::FREE].next = self.entries.len() - 1;
        }
        let index = self.entries[Self::FREE].next;
        self.entries[index].value = value;
        self.unlink(index);
        self.link_after(index, Self::OCCUPIED);

        index
    }

    fn remove(&mut self, index: usize) -> T {
        self.unlink(index);
        self.link_after(index, Self::FREE);

        self.entries[index].value.take().expect("invalid index")
    }

    fn back(&self) -> usize {
        self.entries[Self::OCCUPIED].prev
    }

    fn pop_back(&mut self) -> T {
        let index = self.back();

        self.remove(index)
    }

    fn get(&self, index: usize) -> &T {
        self.entries[index].value.as_ref().expect("invalid index")
    }

    fn get_mut(&mut self, index: usize) -> &mut T {
        self.entries[index].value.as_mut().expect("invalid index")
    }

    fn set(&mut self, index: usize, value: T) -> Option<T> {
        std::mem::replace(&mut self.entries[index].value, Some(value))
    }

    fn clear(&mut self) {
        self.entries.clear();
        self.entries.push(ListEntry::<T> {
            value: None,
            next: 0,
            prev: 0,
        });
        self.entries.push(ListEntry::<T> {
            value: None,
            next: 1,
            prev: 1,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::collections::hash_map::RandomState;
    #[test]
    fn get_after_insert() {
        let mut cache = LruCache::with_hasher(RandomState::default());

        assert_eq!(cache.get(&1), None);

        cache.insert(1, 2);

        assert_eq!(cache.get(&1), Some(&2));
    }
    #[test]
    fn evict_lru() {
        let mut cache = LruCache::with_hasher(RandomState::default());

        cache.insert(1, 2);
        cache.insert(2, 3);
        cache.insert(3, 4);
        cache.insert(4, 5);
        cache.insert(2, 5);
        cache.get(&1);

        assert_eq!(cache.evict_lru(), Some((3, 4)));
        assert_eq!(cache.evict_lru(), Some((4, 5)));
        assert_eq!(cache.evict_lru(), Some((2, 5)));
        assert_eq!(cache.evict_lru(), Some((1, 2)));

        assert!(cache.is_empty());
    }
}
