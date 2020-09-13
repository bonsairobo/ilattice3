use crate::{
    bounding_extent, copy_extent,
    lattice_map::LatticeMapKeyValIterator,
    prelude::*,
    vec_lattice_map::{FastCompressedVecLatticeMap, FastLZ4},
    Extent, Point, VecLatticeMap, YLevelsIndexer,
};

use compressible_map::{
    Compressible, CompressibleMap, Decompressible, LocalCache, MaybeCompressed,
};
use serde::{Deserialize, Serialize};

/// Stores a partial function on ZxZxZ in same-sized chunks using a `CompressibleMap`. Optionally,
/// you can also store some metadata of type `M` for each chunk.
pub struct ChunkedLatticeMap<T, M = (), I = YLevelsIndexer>
where
    I: Indexer,
    M: Clone,
{
    chunk_size: Point,
    pub chunks: CompressibleFnvMap<Point, Chunk<T, M, I>, FastLZ4>,
}

type CompressibleFnvMap<K, V, A> = CompressibleMap<K, V, A, fnv::FnvBuildHasher>;

pub type LocalChunkCache<T, M, I> = LocalCache<Point, Chunk<T, M, I>, fnv::FnvBuildHasher>;

/// One piece of the `ChunkedLatticeMap`. Contains both some generic metadata and the data for each
/// point in the chunk extent.
#[derive(Clone, Deserialize, Serialize)]
pub struct Chunk<T, M = (), I = YLevelsIndexer> {
    pub metadata: M,
    pub map: VecLatticeMap<T, I>,
}

impl<T, I> Chunk<T, (), I> {
    /// Constructs a chunk without metadata.
    pub fn with_map(map: VecLatticeMap<T, I>) -> Self {
        Chunk { metadata: (), map }
    }
}

pub struct FastCompressedChunk<T, M, I> {
    metadata: M, // metadata doesn't get compressed, hope it's small!
    compressed_map: FastCompressedVecLatticeMap<T, I>,
}

// PERF: cloning the metadata is unfortunate

impl<T, M, I> Decompressible<FastLZ4> for FastCompressedChunk<T, M, I>
where
    I: Indexer,
    M: Clone,
{
    type Decompressed = Chunk<T, M, I>;

    fn decompress(&self) -> Self::Decompressed {
        Chunk {
            metadata: self.metadata.clone(),
            map: self.compressed_map.decompress(),
        }
    }
}

impl<T, M, I> Compressible<FastLZ4> for Chunk<T, M, I>
where
    I: Indexer,
    M: Clone,
{
    type Compressed = FastCompressedChunk<T, M, I>;

    fn compress(&self, params: FastLZ4) -> Self::Compressed {
        FastCompressedChunk {
            metadata: self.metadata.clone(),
            compressed_map: self.map.compress(params),
        }
    }
}

impl<T, M, I> ChunkedLatticeMap<T, M, I>
where
    I: Indexer,
    M: Clone,
{
    /// Creates an empty map.
    pub fn new(chunk_size: Point) -> Self {
        assert!(chunk_size > [0, 0, 0].into());

        Self {
            chunk_size,
            // TODO: don't hardcore the compression level
            chunks: CompressibleFnvMap::new(FastLZ4 { level: 10 }),
        }
    }

    pub fn chunk_size(&self) -> &Point {
        &self.chunk_size
    }

    /// Returns the key of the chunk that contains `point`.
    pub fn chunk_key(&self, point: &Point) -> Point {
        *point / self.chunk_size
    }

    /// Returns the extent whose points are exactly the set of chunks (keys) that contain any of the
    /// points in `extent`. For example, if the chunk size is 2x2x2, then the extent from (0, 0, 0)
    /// to (8, 8, 8) would return the chunk key extent from (0, 0, 0) to (4, 4, 4).
    fn key_extent(&self, extent: &Extent) -> Extent {
        let key_min = self.chunk_key(&extent.get_minimum());
        let key_max = self.chunk_key(&extent.get_world_max());

        Extent::from_min_and_world_max(key_min, key_max)
    }

    pub fn extent_for_chunk_key(&self, key: &Point) -> Extent {
        extent_for_chunk_key(&self.chunk_size, key)
    }

    /// Returns the chunk at `key` if it exists.
    pub fn get_chunk<'a>(
        &'a self,
        key: Point,
        local_cache: &'a LocalChunkCache<T, M, I>,
    ) -> Option<&Chunk<T, M, I>> {
        self.chunks.get_const(key, local_cache)
    }

    /// Returns the mutable chunk at `key` if it exists.
    pub fn get_mut_chunk(&mut self, key: Point) -> Option<&mut Chunk<T, M, I>> {
        self.chunks.get_mut(key)
    }

    /// Returns the chunk containing `point` if it exists.
    pub fn get_chunk_containing_point<'a>(
        &'a mut self,
        point: &Point,
        local_cache: &'a LocalChunkCache<T, M, I>,
    ) -> Option<(Point, &Chunk<T, M, I>)> {
        let chunk_key = self.chunk_key(point);

        self.get_chunk(chunk_key, local_cache)
            .map(|c| (chunk_key, c))
    }

    /// Returns the mutable chunk containing `point` if it exists.
    pub fn get_mut_chunk_containing_point(
        &mut self,
        point: &Point,
    ) -> Option<(Point, &mut Chunk<T, M, I>)> {
        let chunk_key = self.chunk_key(point);

        self.get_mut_chunk(chunk_key).map(|c| (chunk_key, c))
    }

    /// Returns an iterator over all points and corresponding values in the given extent. If chunks
    /// are missing from the extent, then their points will not be yielded.
    pub fn iter_point_values<'a>(
        &'a self,
        extent: Extent,
        local_cache: &'a LocalChunkCache<T, M, I>,
    ) -> impl Iterator<Item = (Point, &T)> {
        let key_extent = self.key_extent(&extent);

        let mut chunks = Vec::new();
        for key in key_extent {
            if let Some(chunk) = self.chunks.get_const(key, local_cache) {
                chunks.push(&chunk.map);
            }
        }

        chunks.into_iter().flat_map(move |chunk| {
            let extent_iter = chunk.get_extent().intersection(&extent).into_iter();

            LatticeMapKeyValIterator::new(chunk, extent_iter)
        })
    }

    /// An iterator over the chunk keys (sparse points in the space of chunk coordinates).
    pub fn chunk_keys(&self) -> impl Iterator<Item = &Point> {
        self.chunks.keys()
    }

    /// Returns the smallest extent that bounds all points in this map.
    pub fn bounding_extent(&self) -> Extent {
        // PERF: There's probably a way to avoid iterating over all of the chunks.
        assert!(!self.chunks.is_empty());
        let extrema_iter = self
            .chunks
            .iter_maybe_compressed()
            .flat_map(|chunk| match chunk.1 {
                MaybeCompressed::Decompressed(c) => {
                    let extent = c.map.get_extent();

                    vec![extent.get_minimum(), extent.get_world_max()].into_iter()
                }
                MaybeCompressed::Compressed(c) => {
                    let extent = c.compressed_map.get_extent();

                    vec![extent.get_minimum(), extent.get_world_max()].into_iter()
                }
            });

        bounding_extent(extrema_iter)
    }
}

/// Returns the extent of the chunk at `key`.
pub fn extent_for_chunk_key(chunk_size: &Point, key: &Point) -> Extent {
    let min = *key * *chunk_size;
    let local_sup = *chunk_size;

    Extent::from_min_and_local_supremum(min, local_sup)
}

impl<T: Clone, M, I> ChunkedLatticeMap<T, M, I>
where
    M: Clone,
    I: Indexer,
{
    /// Get mutable data for point `p`. If `p` does not exist, calls `fill_empty_chunk` to fill
    /// that entry first.
    pub fn get_mut_or_create(
        &mut self,
        p: &Point,
        fill_empty_chunk: impl Fn(&Point, &Extent) -> Chunk<T, M, I>,
    ) -> (Point, &mut T) {
        let key = self.chunk_key(p);
        let chunk_size = self.chunk_size;

        (
            key,
            self.chunks
                .get_or_insert_with(key, || {
                    fill_empty_chunk(&key, &extent_for_chunk_key(&chunk_size, &key))
                })
                .map
                .get_world_ref_mut(p),
        )
    }

    /// Fills `extent` with value `T`. If `extent` overlaps a chunk that doesn't exist yet, then
    /// `fill_empty_chunk` will fill the chunk first.
    pub fn fill_extent(
        &mut self,
        extent: &Extent,
        val: T,
        fill_empty_chunk: impl Fn(&Point, &Extent) -> Chunk<T, M, I>,
    ) {
        let fill_lat = VecLatticeMap::<_, YLevelsIndexer>::fill(*extent, val);
        self.copy_map_into_chunks(&fill_lat, fill_empty_chunk);
    }

    /// Copies the values at all points in `map` into the same points of `self`. For any points
    /// inside of chunks that don't exist yet, first fill the chunk with `fill_empty_chunk`.
    pub fn copy_map_into_chunks<V>(
        &mut self,
        map: &V,
        fill_empty_chunk: impl Fn(&Point, &Extent) -> Chunk<T, M, I>,
    ) where
        V: GetWorld<Data = T> + GetExtent,
    {
        for key in self.key_extent(&map.get_extent()) {
            let chunk_extent = self.extent_for_chunk_key(&key);
            let chunk = self
                .chunks
                .get_or_insert_with(key, || fill_empty_chunk(&key, &chunk_extent));
            copy_extent(
                map,
                &mut chunk.map,
                &chunk_extent.intersection(&map.get_extent()),
            );
        }
    }
}

impl<T, M, I> ChunkedLatticeMap<T, M, I>
where
    T: Clone,
    M: Clone,
    I: Indexer,
{
    /// Fills `extent` with value `T`. If `extent` overlaps a chunk that doesn't exist yet, then
    /// `default_voxel` will fill the chunk first.
    pub fn fill_extent_or_default(
        &mut self,
        extent: &Extent,
        val: T,
        default_metadata: M,
        default_voxel: T,
    ) {
        let fill_lat = VecLatticeMap::<_, YLevelsIndexer>::fill(*extent, val);
        self.copy_map_into_chunks(&fill_lat, |_, extent| Chunk {
            metadata: default_metadata.clone(),
            map: VecLatticeMap::fill(*extent, default_voxel.clone()),
        });
    }

    /// Sets point `p` to value `T`. If `p` is in a chunk that doesn't exist yet, then
    /// `default_voxel` will fill the chunk first.
    pub fn get_mut_or_default(
        &mut self,
        p: &Point,
        default_metadata: M,
        default_voxel: T,
    ) -> (Point, &mut T) {
        self.get_mut_or_create(p, |_, extent| Chunk {
            metadata: default_metadata.clone(),
            map: VecLatticeMap::fill(*extent, default_voxel.clone()),
        })
    }
}

impl<T, M, I> ChunkedLatticeMap<T, M, I>
where
    T: Clone + Default,
    M: Clone,
    I: Indexer,
{
    /// Copies all points from `extent` into a new dense map. Any points in `self` that have values
    /// will take on the value `T::default()`.
    pub fn copy_extent_into_new_map(
        &self,
        extent: Extent,
        local_cache: &LocalChunkCache<T, M, I>,
    ) -> VecLatticeMap<T, I> {
        let mut new_map = VecLatticeMap::fill(extent, T::default());
        for (p, val) in self.iter_point_values(extent, local_cache) {
            *new_map.get_world_ref_mut(&p) = val.clone();
        }

        new_map
    }

    /// Copies all points from `self` into a new dense map. Any points in `self` that have values
    /// will take on the value `T::default()`.
    pub fn copy_into_new_map(&self, local_cache: &LocalChunkCache<T, M, I>) -> VecLatticeMap<T, I> {
        self.copy_extent_into_new_map(self.bounding_extent(), local_cache)
    }

    /// Copies the entire chunk and values at adjacent points into a new dense map. For example,
    /// if the key for chunk with extent from (2, 2, 2) to (4, 4, 4) is chose, then the extent from
    /// (1, 1, 1) to (5, 5, 5) will be copied.
    pub fn get_chunk_and_boundary(
        &self,
        chunk_key: &Point,
        local_cache: &LocalChunkCache<T, M, I>,
    ) -> VecLatticeMap<T, I> {
        let extent = self.extent_for_chunk_key(chunk_key).padded(1);

        self.copy_extent_into_new_map(extent, local_cache)
    }
}

impl<T, M, I> MaybeGetWorldRefMut for ChunkedLatticeMap<T, M, I>
where
    M: Clone,
    I: Indexer,
{
    type Data = T;

    fn maybe_get_world_ref_mut(&mut self, p: &Point) -> Option<&mut T> {
        self.get_mut_chunk_containing_point(p)
            .map(|(_key, chunk)| chunk.map.get_world_ref_mut(p))
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
    use crate::test_util::assert_elements_eq;

    #[test]
    fn test_chunked_map_iterator() {
        // All points less than [0, 0, 0] in the query extent will have value 1.
        let submap1 = VecLatticeMap::<_, YLevelsIndexer>::fill(
            Extent::from_center_and_radius([-7, -7, -7].into(), 7),
            1,
        );
        // Only [1, 1, 1] with have value 2 in the query extent.
        let submap2 = VecLatticeMap::<_, YLevelsIndexer>::fill(
            Extent::from_center_and_radius([8, 8, 8].into(), 7),
            2,
        );
        // All other point values will not change (stay 0).

        let mut chunked_map: ChunkedLatticeMap<u32> = ChunkedLatticeMap::new([4, 4, 4].into());
        let fill_empty_chunk =
            |_key: &Point, extent: &Extent| Chunk::with_map(VecLatticeMap::fill(*extent, 0));
        chunked_map.copy_map_into_chunks(&submap1, &fill_empty_chunk);
        chunked_map.copy_map_into_chunks(&submap2, &fill_empty_chunk);

        let query_extent = Extent::from_center_and_radius([0, 0, 0].into(), 1);
        let local_cache = LocalChunkCache::new();
        let points: Vec<_> = chunked_map
            .iter_point_values(query_extent, &local_cache)
            .map(|(p, i)| (p, *i))
            .collect();

        assert_elements_eq(
            &points,
            &vec![
                ([-1, -1, -1].into(), 1),
                ([-1, -1, 0].into(), 1),
                ([-1, -1, 1].into(), 0),
                ([-1, 0, -1].into(), 1),
                ([-1, 0, 0].into(), 1),
                ([-1, 0, 1].into(), 0),
                ([-1, 1, -1].into(), 0),
                ([-1, 1, 0].into(), 0),
                ([-1, 1, 1].into(), 0),
                ([0, -1, -1].into(), 1),
                ([0, -1, 0].into(), 1),
                ([0, -1, 1].into(), 0),
                ([0, 0, -1].into(), 1),
                ([0, 0, 0].into(), 1),
                ([0, 0, 1].into(), 0),
                ([0, 1, -1].into(), 0),
                ([0, 1, 0].into(), 0),
                ([0, 1, 1].into(), 0),
                ([1, -1, -1].into(), 0),
                ([1, -1, 0].into(), 0),
                ([1, -1, 1].into(), 0),
                ([1, 0, -1].into(), 0),
                ([1, 0, 0].into(), 0),
                ([1, 0, 1].into(), 0),
                ([1, 1, -1].into(), 0),
                ([1, 1, 0].into(), 0),
                ([1, 1, 1].into(), 2),
            ],
        );
    }
}
