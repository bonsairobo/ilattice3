use crate::{
    prelude::*, ChunkedLatticeMap, Extent, HasIndexer, Indexer, Point, VecLatticeMap,
    YLevelsIndexer,
};
use serde::{Deserialize, Serialize};

/// An owned `VecLatticeMap` with a borrowed palette. The `VecLatticeMap` does not need to be a
/// chunk.
pub struct LatticeVoxels<'a, T, P, I> {
    infos: &'a Vec<T>,
    pub map: VecLatticeMap<P, I>,
}

impl<'a, T, P, I> LatticeVoxels<'a, T, P, I>
where
    P: Into<usize>,
{
    pub fn get_pointed_voxel_info(&'a self, ptr: P) -> &'a T {
        &self.infos[ptr.into()]
    }
}

impl<'a, T, P, I> GetExtent for LatticeVoxels<'a, T, P, I>
where
    P: Into<usize>,
{
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

// TODO: impl GetLinear* for LatticeVoxels

impl<'a, T, P, I> GetLinearRef for LatticeVoxels<'a, T, P, I>
where
    P: Clone + Into<usize>,
    I: Indexer,
{
    type Data = T;

    fn get_linear_ref(&self, i: usize) -> &T {
        self.get_pointed_voxel_info(self.map.get_linear(i))
    }
}

/// A borrowed chunk and palette.
pub struct ChunkVoxelsRef<'a, T, P, I> {
    infos: &'a Vec<T>,
    pub map: &'a VecLatticeMap<P, I>,
}

impl<'a, T, P, I> ChunkVoxelsRef<'a, T, P, I>
where
    P: Into<usize>,
{
    pub fn get_pointed_voxel_info(&self, ptr: P) -> &T {
        &self.infos[ptr.into()]
    }
}

impl<'a, T, P, I> GetExtent for ChunkVoxelsRef<'a, T, P, I> {
    fn get_extent(&self) -> &Extent {
        self.map.get_extent()
    }
}

// TODO: impl HasIndexer and GetLinear* for ChunkVoxels*

impl<'a, T, P, I> GetLocalRef for ChunkVoxelsRef<'a, T, P, I>
where
    P: Clone + Into<usize>,
    I: Indexer,
{
    type Data = T;

    fn get_local_ref(&self, p: &Point) -> &T {
        self.get_pointed_voxel_info(self.map.get_local(p))
    }
}

/// A mutably borrowed chunk and palette.
pub struct ChunkVoxelsRefMut<'a, T, P, I> {
    infos: &'a mut Vec<T>,
    map: &'a mut VecLatticeMap<P, I>,
}

impl<'a, T, P, I> ChunkVoxelsRefMut<'a, T, P, I>
where
    P: Into<usize>,
{
    pub fn get_pointed_voxel_info_mut(&mut self, ptr: P) -> &mut T {
        &mut self.infos[ptr.into()]
    }
}

impl<'a, T, P, I> GetExtent for ChunkVoxelsRefMut<'a, T, P, I> {
    fn get_extent(&self) -> &Extent {
        self.map.get_extent()
    }
}

impl<'a, T, P, I> GetLocalRefMut for ChunkVoxelsRefMut<'a, T, P, I>
where
    P: Clone + Into<usize>,
    I: Indexer,
{
    type Data = T;

    fn get_local_ref_mut(&mut self, p: &Point) -> &mut T {
        self.get_pointed_voxel_info_mut(self.map.get_local(p))
    }
}

/// A `ChunkedLatticeMap` that stores voxel data efficiently using palette compression.
#[derive(Deserialize, Serialize)]
pub struct ChunkedPaletteLatticeMap<T, P, I = YLevelsIndexer> {
    /// The palette of voxels that can be used in the map.
    pub infos: Vec<T>,
    /// Which voxels are used at specific points of the lattice.
    pub map: ChunkedLatticeMap<P, I>,
}

impl<T, P, I> ChunkedPaletteLatticeMap<T, P, I> {
    pub fn get_chunk(&self, chunk_key: &Point) -> Option<ChunkVoxelsRef<T, P, I>> {
        let ChunkedPaletteLatticeMap { map, infos } = self;

        map.get_chunk(chunk_key).map(|chunk_map| ChunkVoxelsRef {
            map: chunk_map,
            infos,
        })
    }

    pub fn get_chunk_mut(&mut self, chunk_key: &Point) -> Option<ChunkVoxelsRefMut<T, P, I>> {
        let ChunkedPaletteLatticeMap { map, infos } = self;

        map.get_mut_chunk(chunk_key)
            .map(move |chunk_map| ChunkVoxelsRefMut {
                map: chunk_map,
                infos,
            })
    }

    pub fn iter_chunks_ref(&self) -> impl Iterator<Item = (&Point, ChunkVoxelsRef<T, P, I>)> {
        self.map.iter_chunks().map(move |(chunk_key, chunk_map)| {
            (
                chunk_key,
                ChunkVoxelsRef {
                    map: chunk_map,
                    infos: &self.infos,
                },
            )
        })
    }
}

impl<T, P, I> ChunkedPaletteLatticeMap<T, P, I>
where
    P: Clone + Default,
    I: Indexer,
{
    pub fn get_chunk_and_boundary(&self, chunk_key: &Point) -> LatticeVoxels<T, P, I> {
        let ChunkedPaletteLatticeMap { map, infos } = self;

        let chunk_and_boundary = map.get_chunk_and_boundary(chunk_key);

        LatticeVoxels {
            infos,
            map: chunk_and_boundary,
        }
    }

    pub fn iter_chunks_with_boundary(
        &self,
    ) -> impl Iterator<Item = (&Point, LatticeVoxels<T, P, I>)> {
        self.map
            .chunk_keys()
            .map(move |chunk_key| (chunk_key, self.get_chunk_and_boundary(chunk_key)))
    }
}

impl<T, P> MaybeGetWorldRef for ChunkedPaletteLatticeMap<T, P>
where
    P: Clone + Into<usize>,
{
    type Data = T;

    fn maybe_get_world_ref(&self, p: &Point) -> Option<&T> {
        self.map
            .maybe_get_world(p)
            .map(|ptr| &self.infos[ptr.into()])
    }
}

impl<T, P> MaybeGetWorldRefMut for ChunkedPaletteLatticeMap<T, P>
where
    P: Clone + Into<usize>,
{
    type Data = T;

    fn maybe_get_world_ref_mut(&mut self, p: &Point) -> Option<&mut T> {
        self.map
            .maybe_get_world(p)
            .map(move |ptr| &mut self.infos[ptr.into()])
    }
}
