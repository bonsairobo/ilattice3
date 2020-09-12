use crate::{
    lru_cache::{EntryState, LruCache},
    Compressible, Decompressible,
};

use fnv::FnvBuildHasher;
use std::collections::HashMap;
use std::hash::{BuildHasher, Hash};

// TODO: we could probably reuse compressed values if the decompressed value doesn't get modified

/// A hash map that allows compressing the least recently used values. Mostly useful when you need
/// to store a lot of large values in memory.
///
/// Call the `compress_lru` method to compress the least recently used value. The most recently used
/// values will stay uncompressed in a cache.
///
/// Any **mutable** access (`&mut self`) that misses the cache will decompress and cache the value
/// inline. You can call `get` to prefetch into the cache and avoid extra latency on further
/// accesses.
///
/// Any **immutable** access (`&self`, e.g. from multiple threads), like `get_const`, cannot update
/// the cache. Instead, it will record accesses and store decompressed values in a
/// `ThreadLocalCache` that can be used later to update the cache with `flush_thread_local_cache`.
pub struct CompressibleMap<K, V, H, A>
where
    V: Compressible<A>,
{
    cache: LruCache<K, V, H>,
    compressed: HashMap<K, <V as Compressible<A>>::Compressed, H>,
    compression_params: A,
}

#[derive(Default)]
pub struct ThreadLocalCache<K, V, H> {
    accesses: HashMap<K, Option<V>, H>,
}

pub type CompressibleFnvMap<K, V, A> = CompressibleMap<K, V, FnvBuildHasher, A>;

impl<K, V, Vc, A> CompressibleFnvMap<K, V, A>
where
    K: Clone + Eq + Hash,
    V: Compressible<A, Compressed = Vc>,
    Vc: Decompressible<A, Decompressed = V>,
    A: Clone,
{
    pub fn new_fnv(compression_params: A) -> Self {
        CompressibleMap::new(compression_params)
    }
}

impl<K, V, Vc, H, A> CompressibleMap<K, V, H, A>
where
    K: Clone + Eq + Hash,
    V: Compressible<A, Compressed = Vc>,
    Vc: Decompressible<A, Decompressed = V>,
    H: BuildHasher + Default,
    A: Clone,
{
    pub fn new(compression_params: A) -> Self {
        Self {
            cache: LruCache::default(),
            compressed: HashMap::default(),
            compression_params,
        }
    }

    /// Insert a new value and drop the old one.
    pub fn insert(&mut self, key: K, value: V) {
        self.cache.insert(key.clone(), value);

        // PERF: this might not be necessary, but we need to confirm that the compressed value won't
        // pop up again somewhere and cause inconsistencies
        self.compressed.remove(&key);
    }

    /// Insert a new value and return the old one if it exists, which requires decompressing it.
    pub fn replace(&mut self, key: K, value: V) -> Option<V> {
        self.cache
            .insert(key.clone(), value)
            .map(|old_cache_entry| match old_cache_entry {
                EntryState::Cached(v) => v,
                EntryState::Evicted => {
                    let compressed_value = self.compressed.remove(&key).unwrap();

                    compressed_value.decompress()
                }
            })
    }

    pub fn compress_lru(&mut self) {
        if let Some((lru_key, lru_value)) = self.cache.evict_lru() {
            self.compressed
                .insert(lru_key, lru_value.compress(self.compression_params.clone()));
        }
    }

    pub fn get_mut(&mut self, key: K) -> Option<&mut V> {
        let CompressibleMap {
            cache, compressed, ..
        } = self;

        cache.get_or_repopulate_with(key.clone(), || {
            compressed.remove(&key).map(|v| v.decompress()).unwrap()
        })
    }

    pub fn get(&mut self, key: K) -> Option<&V> {
        // Hopefully downgrading the reference is a NOOP.
        self.get_mut(key).map(|v| &*v)
    }

    /// Used for thread-safe access. The cache will not be updated, but accesses will be recorded in
    /// the provided `ThreadLocalCache`. Call `flush_thread_local_cache` to update the cache with
    /// the record.
    pub fn get_const<'a>(
        &'a self,
        key: K,
        local_cache: &'a mut ThreadLocalCache<K, V, H>,
    ) -> Option<&'a V> {
        self.cache.get_const(&key).map(move |entry| {
            match entry {
                EntryState::Cached(v) => {
                    // For the sake of updating LRU order when we flush this local cache.
                    local_cache.accesses.insert(key, None);

                    v
                }
                EntryState::Evicted => {
                    // Check the thread-local record before trying to decompress.
                    local_cache
                        .accesses
                        .entry(key.clone())
                        .or_insert_with(move || {
                            // OK, we have to decompress.
                            Some(self.compressed.get(&key).unwrap().decompress())
                        })
                        .as_ref()
                        .unwrap()
                }
            }
        })
    }

    /// Updates the cache and it's approximate LRU order.
    pub fn flush_thread_local_cache(&mut self, record: ThreadLocalCache<K, V, H>) {
        for (key, value) in record.accesses.into_iter() {
            if let Some(value) = value {
                self.insert(key, value);
            } else {
                // None value means there was a cached access
                self.cache.get(&key);
            }
        }
    }

    pub fn drop(&mut self, key: &K) {
        self.cache.remove(key);
        self.compressed.remove(key);
    }

    /// Removes the value and returns it if it exists, decompressing it first.
    pub fn remove(&mut self, key: &K) -> Option<V> {
        self.cache.remove(key).map(|entry| match entry {
            EntryState::Cached(v) => v,
            EntryState::Evicted => self.compressed.remove(key).unwrap().decompress(),
        })
    }

    pub fn clear(&mut self) {
        self.cache.clear();
        self.compressed.clear();
    }

    pub fn len(&self) -> usize {
        self.len_cached() + self.compressed_len()
    }

    pub fn len_cached(&self) -> usize {
        self.cache.len_cached()
    }

    pub fn compressed_len(&self) -> usize {
        self.compressed.len()
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

    #[derive(Debug, Default, Eq, PartialEq)]
    struct Foo(u32);

    struct FooCompressed(u32);

    impl Compressible<()> for Foo {
        type Compressed = FooCompressed;

        fn compress(&self, _: ()) -> Self::Compressed {
            FooCompressed(self.0 + 1)
        }
    }

    impl Decompressible<()> for FooCompressed {
        type Decompressed = Foo;

        fn decompress(&self) -> Self::Decompressed {
            Foo(self.0 + 1)
        }
    }

    #[test]
    fn get_after_compress() {
        let mut map = CompressibleFnvMap::new_fnv(());

        map.insert(1, Foo(0));

        map.compress_lru();

        assert_eq!(map.len_cached(), 0);
        assert_eq!(map.compressed_len(), 1);

        assert_eq!(Some(&Foo(2)), map.get(1));

        assert_eq!(map.len_cached(), 1);
        assert_eq!(map.compressed_len(), 0);
    }

    #[test]
    fn flush_after_get_const_populates_cache() {
        let mut map = CompressibleFnvMap::new_fnv(());

        map.insert(1, Foo(0));

        map.compress_lru();

        let mut record = ThreadLocalCache::default();
        assert_eq!(Some(&Foo(2)), map.get_const(1, &mut record));

        assert_eq!(map.len_cached(), 0);
        assert_eq!(map.compressed_len(), 1);

        map.flush_thread_local_cache(record);

        assert_eq!(map.len_cached(), 1);
        assert_eq!(map.compressed_len(), 0);

        assert_eq!(Some(&Foo(2)), map.get(1));
    }
}
