//! It is customary for users to store data in two places: the palette and the underlying lattice
//! map. Users may want to implement the Get* traits in multiple ways, and they can use newtypes
//! that wrap these types to do so.

use crate::{
    chunked_lattice_map::SerializableChunkedLatticeMap, prelude::*, vec_lattice_map::FastLZ4,
    ChunkedLatticeMap, Extent, LocalChunkCache, Point, VecLatticeMap, YLevelsIndexer,
};

use compressible_map::BincodeLz4;
use serde::{de::DeserializeOwned, Deserialize, Serialize};

pub trait GetPaletteAddress {
    fn get_palette_address(&self) -> usize;
}

/// A `ChunkedLatticeMap` that stores *shared* voxel data efficiently using palette compression. The
/// palette pointer type `P` is what's stored for each point of the lattice. That pointer must
/// implement `GetPaletteAddress` because it is used to look up a `T` value in the palette.
pub struct PaletteLatticeMap<T, P, M = (), I = YLevelsIndexer>
where
    M: Clone,
    I: Indexer,
{
    /// The palette of voxels that can be used in the map.
    pub palette: Vec<T>,
    /// Which voxels are used at specific points of the lattice.
    pub map: ChunkedLatticeMap<P, M, I>,
}

impl<T, P, M, I> PaletteLatticeMap<T, P, M, I>
where
    M: Clone,
    I: Indexer,
{
    /// Returns a chunk reference along with a palette reference.
    pub fn get_chunk<'a>(
        &'a self,
        chunk_key: Point,
        local_cache: &'a LocalChunkCache<P, M, I>,
    ) -> Option<ChunkVoxelsRef<T, P, I>> {
        let PaletteLatticeMap { map, palette } = self;

        map.get_chunk(chunk_key, local_cache)
            .map(|chunk| ChunkVoxelsRef {
                map: &chunk.map,
                palette,
            })
    }

    /// Returns a chunk reference along with a palette reference.
    pub fn get_chunk_containing_point<'a>(
        &'a self,
        point: Point,
        local_cache: &'a LocalChunkCache<P, M, I>,
    ) -> Option<(Point, ChunkVoxelsRef<T, P, I>)> {
        let chunk_key = self.map.chunk_key(&point);

        self.get_chunk(chunk_key, local_cache)
            .map(|chunk| (chunk_key, chunk))
    }

    /// Returns a mutable chunk reference along with a mutable palette reference.
    pub fn get_chunk_mut(&mut self, chunk_key: Point) -> Option<ChunkVoxelsRefMut<T, P, I>> {
        let PaletteLatticeMap { map, palette } = self;

        map.get_mut_chunk(chunk_key)
            .map(move |chunk| ChunkVoxelsRefMut {
                map: &mut chunk.map,
                palette,
            })
    }
}

impl<T, P, M, I> PaletteLatticeMap<T, P, M, I>
where
    P: Clone + Default,
    M: Clone,
    I: Indexer,
{
    /// Like `ChunkedLatticeMap::copy_extent_into_new_map`, but also returns a palette reference.
    pub fn copy_extent_into_new_map(
        &self,
        extent: Extent,
        local_cache: &LocalChunkCache<P, M, I>,
    ) -> LatticeVoxels<T, P, I> {
        let PaletteLatticeMap { map, palette } = self;

        let voxels = map.copy_extent_into_new_map(extent, local_cache);

        LatticeVoxels {
            palette,
            map: voxels,
        }
    }

    /// Like `ChunkedLatticeMap::get_chunk_and_boundary`, but also returns a palette reference.
    pub fn get_chunk_and_boundary(
        &self,
        chunk_key: &Point,
        local_cache: &LocalChunkCache<P, M, I>,
    ) -> LatticeVoxels<T, P, I> {
        let PaletteLatticeMap { map, palette } = self;

        let chunk_and_boundary = map.get_chunk_and_boundary(chunk_key, local_cache);

        LatticeVoxels {
            palette,
            map: chunk_and_boundary,
        }
    }

    /// Iterate over (chunk key, owned lattice map + palette reference) pairs.
    pub fn iter_chunks_with_boundary<'a>(
        &'a self,
        local_cache: &'a LocalChunkCache<P, M, I>,
    ) -> impl Iterator<Item = (&Point, LatticeVoxels<T, P, I>)> {
        self.map.chunk_keys().map(move |chunk_key| {
            (
                chunk_key,
                self.get_chunk_and_boundary(chunk_key, local_cache),
            )
        })
    }
}

#[derive(Deserialize, Serialize)]
pub struct SerializablePaletteLatticeMap<T, P, M, I> {
    pub compressed_chunks: SerializableChunkedLatticeMap<P, M, I>,
    pub palette: Vec<T>,
}

impl<T, P, M, I> PaletteLatticeMap<T, P, M, I>
where
    T: Clone + DeserializeOwned + Serialize,
    P: DeserializeOwned + Serialize,
    M: Clone + DeserializeOwned + Serialize,
    I: Indexer + DeserializeOwned + Serialize,
{
    pub fn to_serializable(&self, params: BincodeLz4) -> SerializablePaletteLatticeMap<T, P, M, I> {
        SerializablePaletteLatticeMap {
            compressed_chunks: self.map.to_serializable(params),
            palette: self.palette.clone(),
        }
    }

    pub fn from_serializable(
        map: SerializablePaletteLatticeMap<T, P, M, I>,
        params: FastLZ4,
    ) -> Self {
        Self {
            palette: map.palette,
            map: ChunkedLatticeMap::from_serializable(&map.compressed_chunks, params),
        }
    }
}

impl<T, P, M, I> MaybeGetWorldRefMut for PaletteLatticeMap<T, P, M, I>
where
    P: Clone + GetPaletteAddress,
    M: Clone,
    I: Indexer,
{
    type Data = T;

    fn maybe_get_world_ref_mut(&mut self, p: &Point) -> Option<&mut T> {
        self.map
            .maybe_get_world_ref_mut(p)
            .cloned()
            .map(move |ptr| &mut self.palette[ptr.get_palette_address()])
    }
}

/// A thread-local reader of a `PaletteLatticeMap` which stores a cache of chunks that were
/// decompressed after missing the global cache of chunks.
pub struct PaletteLatticeMapReader<'a, T, P, M, I>
where
    M: Clone,
    I: Clone + Indexer,
{
    pub map: &'a PaletteLatticeMap<T, P, M, I>,
    pub local_cache: LocalChunkCache<P, M, I>,
}

impl<'a, T, P, M, I> PaletteLatticeMapReader<'a, T, P, M, I>
where
    M: Clone,
    I: Clone + Indexer,
{
    pub fn new(map: &'a PaletteLatticeMap<T, P, M, I>) -> Self {
        Self {
            map,
            local_cache: LocalChunkCache::new(),
        }
    }
}

impl<'a, T, P, M, I> MaybeGetWorldRef for PaletteLatticeMapReader<'a, T, P, M, I>
where
    P: Clone + GetPaletteAddress,
    M: Clone,
    I: Clone + Indexer,
{
    type Data = T;

    fn maybe_get_world_ref(&self, p: &Point) -> Option<&T> {
        self.map
            .get_chunk_containing_point(*p, &self.local_cache)
            .map(|(_key, chunk)| &chunk.palette[chunk.map.get_world_ref(p).get_palette_address()])
    }
}

/// An owned `VecLatticeMap` with a borrowed palette. The `VecLatticeMap` does not need to be a
/// chunk.
pub struct LatticeVoxels<'a, T, P, I = YLevelsIndexer> {
    pub palette: &'a Vec<T>,
    pub map: VecLatticeMap<P, I>,
}

impl<'a, T, P, I> GetExtent for LatticeVoxels<'a, T, P, I> {
    fn get_extent(&self) -> &Extent {
        self.map.get_extent()
    }
}

impl<'a, T, P, I> HasIndexer for LatticeVoxels<'a, T, P, I>
where
    I: Indexer,
{
    type Indexer = I;
}

impl<'a, T, P, I> GetLinear for LatticeVoxels<'a, T, P, I>
where
    T: Clone,
    P: GetPaletteAddress,
    I: Indexer,
{
    type Data = T;

    fn get_linear(&self, i: usize) -> T {
        self.palette[self.map.get_linear_ref(i).get_palette_address()].clone()
    }
}

impl<'a, T, P, I> GetLinearRef for LatticeVoxels<'a, T, P, I>
where
    P: GetPaletteAddress,
    I: Indexer,
{
    type Data = T;

    fn get_linear_ref(&self, i: usize) -> &T {
        &self.palette[self.map.get_linear_ref(i).get_palette_address()]
    }
}

/// A borrowed chunk and palette.
pub struct ChunkVoxelsRef<'a, T, P, I> {
    pub palette: &'a Vec<T>,
    pub map: &'a VecLatticeMap<P, I>,
}

impl<'a, T, P, I> GetExtent for ChunkVoxelsRef<'a, T, P, I> {
    fn get_extent(&self) -> &Extent {
        self.map.get_extent()
    }
}

impl<'a, T, P, I> HasIndexer for ChunkVoxelsRef<'a, T, P, I>
where
    I: Indexer,
{
    type Indexer = I;
}

impl<'a, T, P, I> GetLinearRef for ChunkVoxelsRef<'a, T, P, I>
where
    P: Clone + GetPaletteAddress,
    I: Indexer,
{
    type Data = T;

    fn get_linear_ref(&self, i: usize) -> &Self::Data {
        &self.palette[self.map.get_linear_ref(i).get_palette_address()]
    }
}

/// A mutably borrowed chunk and palette.
pub struct ChunkVoxelsRefMut<'a, T, P, I> {
    pub palette: &'a mut Vec<T>,
    pub map: &'a mut VecLatticeMap<P, I>,
}

impl<'a, T, P, I> GetExtent for ChunkVoxelsRefMut<'a, T, P, I> {
    fn get_extent(&self) -> &Extent {
        self.map.get_extent()
    }
}

impl<'a, T, P, I> HasIndexer for ChunkVoxelsRefMut<'a, T, P, I>
where
    I: Indexer,
{
    type Indexer = I;
}

impl<'a, T, P, I> GetLinearRefMut for ChunkVoxelsRefMut<'a, T, P, I>
where
    P: Clone + GetPaletteAddress,
    I: Indexer,
{
    type Data = T;

    fn get_linear_ref_mut(&mut self, i: usize) -> &mut Self::Data {
        &mut self.palette[self.map.get_linear_ref(i).get_palette_address()]
    }
}
