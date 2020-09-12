use core::hash::{BuildHasher, Hash};
use std::collections::{hash_map, HashMap};

/// A cache that tracks the Least Recently Used element for next eviction. Here, "used" means read
/// or written.
///
/// Eviction does not happen inline; the user must explicitly call `evict_lru` to evict the LRU
/// element. Thus the cache may grow unbounded unless evictions or explicit removals occur.
///
/// Note that eviction and removal are not treated the same. The cache remembers elements that have
/// been evicted but not removed. This is useful when users need to store evicted data in a separate
/// structure, since if they look up a key and get `Some(EntryState::Evicted)`, they know that the
/// data exists somewhere else. If they get `None`, then they don't have to look elsewhere; the data
/// simply doesn't exist anywhere.
#[derive(Clone, Debug)]
pub struct LruCache<K, V, H> {
    store: HashMap<K, EntryState<usize>, H>,
    order: LruList<(K, V)>,
    num_evicted: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EntryState<V> {
    Cached(V),
    Evicted,
}

impl<V> EntryState<V> {
    pub fn map<T>(self, f: impl FnOnce(V) -> T) -> EntryState<T> {
        match self {
            EntryState::Cached(v) => EntryState::Cached(f(v)),
            EntryState::Evicted => EntryState::Evicted,
        }
    }

    pub fn some_if_cached(self) -> Option<V> {
        match self {
            EntryState::Cached(v) => Some(v),
            EntryState::Evicted => None,
        }
    }
}

impl<K, V, H> Default for LruCache<K, V, H>
where
    H: Default,
    K: Hash + Eq,
{
    fn default() -> Self {
        Self::with_hasher(Default::default())
    }
}

impl<K, V, H> LruCache<K, V, H>
where
    K: Hash + Eq,
{
    pub fn with_hasher(hasher_builder: H) -> LruCache<K, V, H> {
        LruCache {
            store: HashMap::with_hasher(hasher_builder),
            order: LruList::new(),
            num_evicted: 0,
        }
    }
}

impl<K, V, H> LruCache<K, V, H>
where
    K: Hash + Eq + Clone,
    H: BuildHasher,
{
    pub fn get_mut(&mut self, key: &K) -> Option<EntryState<&mut V>> {
        let Self { store, order, .. } = self;
        store.get(key).map(move |&entry| {
            entry.map(move |index| {
                order.move_to_front(index);

                &mut order.get_mut(index).1
            })
        })
    }

    pub fn get(&mut self, key: &K) -> Option<EntryState<&V>> {
        // Hopefully downgrading the reference is a NOOP.
        self.get_mut(key).map(|e| e.map(|v| &*v))
    }

    /// Allows us to get a const reference without having `&mut self`. WARNING: This will not update
    /// the LRU order!
    pub fn get_const(&self, key: &K) -> Option<EntryState<&V>> {
        self.store
            .get(key)
            .cloned()
            .map(|entry| entry.map(|index| &self.order.get(index).1))
    }

    /// Inserts a new `val` for `key`, returning the old entry if it exists.
    pub fn insert(&mut self, key: K, val: V) -> Option<EntryState<V>> {
        let Self { store, order, .. } = self;
        match store.entry(key.clone()) {
            hash_map::Entry::Occupied(mut occupied) => match *occupied.get() {
                EntryState::Cached(index) => {
                    order.move_to_front(index);

                    order
                        .set(index, (key, val))
                        .map(|(_, v)| EntryState::Cached(v))
                }
                EntryState::Evicted => {
                    let new_index = order.push_front(Some((key, val)));
                    occupied.insert(EntryState::Cached(new_index));
                    self.num_evicted -= 1;

                    Some(EntryState::Evicted)
                }
            },
            hash_map::Entry::Vacant(vacant) => {
                vacant.insert(EntryState::Cached(order.push_front(Some((key, val)))));

                None
            }
        }
    }

    /// Tries to get the value for `key`, returning it if it exists. If the entry state is evicted,
    /// calls `f` to repopulate the entry. Otherwise, returns `None`.
    pub fn get_or_repopulate_with(&mut self, key: K, f: impl FnOnce() -> V) -> Option<&mut V> {
        let Self { store, order, .. } = self;
        match store.entry(key.clone()) {
            hash_map::Entry::Occupied(mut occupied) => match *occupied.get() {
                EntryState::Cached(index) => {
                    order.move_to_front(index);

                    Some(&mut order.get_mut(index).1)
                }
                EntryState::Evicted => {
                    let new_index = order.push_front(Some((key, f())));
                    occupied.insert(EntryState::Cached(new_index));
                    self.num_evicted -= 1;

                    Some(&mut order.get_mut(new_index).1)
                }
            },
            hash_map::Entry::Vacant(_) => None,
        }
    }

    // Removes any trace of `key`, such that further accesses will return `None` until a new value
    // is inserted.
    pub fn remove(&mut self, key: &K) -> Option<EntryState<V>> {
        self.store.remove(key).map(|entry| match entry {
            EntryState::Cached(index) => EntryState::Cached(self.order.remove(index).1),
            EntryState::Evicted => {
                self.num_evicted -= 1;

                EntryState::Evicted
            }
        })
    }

    /// Evicts the least-recently used value. This will leave a sentinel behind so that further
    /// accesses will return `Some(EntryState::Evicted)` until the key is removed or a new entry is
    /// inserted.
    pub fn evict_lru(&mut self) -> Option<(K, V)> {
        if self.len_cached() == 0 {
            return None;
        }

        let (key, value) = self.order.pop_back();
        *self.store.get_mut(&key).unwrap() = EntryState::Evicted;
        self.num_evicted += 1;

        Some((key, value))
    }

    /// Removes the least-recently used value, leaving no trace.
    pub fn remove_lru(&mut self) -> Option<(K, V)> {
        if self.len_cached() == 0 {
            return None;
        }

        let (key, value) = self.order.pop_back();
        self.store.remove(&key).unwrap();

        Some((key, value))
    }

    pub fn clear(&mut self) {
        self.store.clear();
        self.order.clear();
        self.num_evicted = 0;
    }

    pub fn len_cached(&self) -> usize {
        self.len_tracked() - self.len_evicted()
    }

    pub fn len_evicted(&self) -> usize {
        self.num_evicted
    }

    pub fn len_tracked(&self) -> usize {
        self.store.len()
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

// ████████╗███████╗███████╗████████╗███████╗
// ╚══██╔══╝██╔════╝██╔════╝╚══██╔══╝██╔════╝
//    ██║   █████╗  ███████╗   ██║   ███████╗
//    ██║   ██╔══╝  ╚════██║   ██║   ╚════██║
//    ██║   ███████╗███████║   ██║   ███████║
//    ╚═╝   ╚══════╝╚══════╝   ╚═╝   ╚══════╝

#[cfg(test)]
mod tests {
    use super::*;

    use std::collections::hash_map::RandomState;

    #[test]
    fn get_after_insert_and_evict_and_remove() {
        let mut cache = LruCache::with_hasher(RandomState::default());
        assert_eq!(cache.get(&1), None);

        cache.insert(1, 2);
        assert_eq!(cache.get(&1), Some(EntryState::Cached(&2)));
        assert_eq!(cache.len_cached(), 1);

        cache.evict_lru();
        assert_eq!(cache.get(&1), Some(EntryState::Evicted));
        assert_eq!(cache.len_evicted(), 1);

        cache.remove(&1);
        assert_eq!(cache.get(&1), None);
        assert_eq!(cache.len_evicted(), 0);
    }

    #[test]
    fn get_after_insert_and_remove() {
        let mut cache = LruCache::with_hasher(RandomState::default());

        cache.insert(1, 2);

        cache.remove(&1);
        assert_eq!(cache.get(&1), None);
        assert_eq!(cache.len_evicted(), 0);
    }

    #[test]
    fn repopulate_after_evict() {
        let mut cache = LruCache::with_hasher(RandomState::default());
        assert_eq!(cache.get(&1), None);

        cache.insert(1, 2);
        cache.evict_lru();

        assert_eq!(cache.get_or_repopulate_with(1, || 3), Some(&mut 3));
        assert_eq!(cache.len_evicted(), 0);
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

        assert!(cache.len_cached() == 0);
        assert!(cache.len_evicted() == 4);
    }

    #[test]
    fn get_const_does_not_affect_lru_order() {
        let mut cache = LruCache::with_hasher(RandomState::default());

        cache.insert(1, 2);
        cache.insert(2, 3);
        cache.get_const(&1);

        assert_eq!(cache.evict_lru(), Some((1, 2)));
        assert_eq!(cache.evict_lru(), Some((2, 3)));
    }
}
