use crate::{prelude::*, ChunkedLatticeMap, Extent, Point, VecLatticeMap};
use serde::{Deserialize, Serialize};

/// An owned `VecLatticeMap` with a borrowed palette. The `VecLatticeMap` does not need to be a
/// chunk.
pub struct LatticeVoxels<'a, T, P> {
    infos: &'a Vec<T>,
    pub map: VecLatticeMap<P>,
}

impl<'a, T, P> LatticeVoxels<'a, T, P>
where
    P: Into<usize>,
{
    pub fn get_pointed_voxel_info(&'a self, ptr: P) -> &'a T {
        &self.infos[ptr.into()]
    }
}

impl<'a, T, P> GetExtent for LatticeVoxels<'a, T, P>
where
    P: Into<usize>,
{
    fn get_extent(&self) -> &Extent {
        self.map.get_extent()
    }
}

impl<'a, T, P> GetLocalRef<T> for LatticeVoxels<'a, T, P>
where
    P: Clone + Into<usize>,
{
    fn get_local_ref(&self, p: &Point) -> &T {
        self.get_pointed_voxel_info(self.map.get_local(p))
    }
}

/// A borrowed chunk and palette.
pub struct ChunkVoxelsRef<'a, T, P> {
    infos: &'a Vec<T>,
    pub map: &'a VecLatticeMap<P>,
}

impl<'a, T, P> ChunkVoxelsRef<'a, T, P>
where
    P: Into<usize>,
{
    pub fn get_pointed_voxel_info(&self, ptr: P) -> &T {
        &self.infos[ptr.into()]
    }
}

impl<'a, T, P> GetExtent for ChunkVoxelsRef<'a, T, P> {
    fn get_extent(&self) -> &Extent {
        self.map.get_extent()
    }
}

impl<'a, T, P> GetLocalRef<T> for ChunkVoxelsRef<'a, T, P>
where
    P: Clone + Into<usize>,
{
    fn get_local_ref(&self, p: &Point) -> &T {
        self.get_pointed_voxel_info(self.map.get_local(p))
    }
}

/// A mutably borrowed chunk and palette.
pub struct ChunkVoxelsRefMut<'a, T, P> {
    infos: &'a mut Vec<T>,
    map: &'a mut VecLatticeMap<P>,
}

impl<'a, T, P> ChunkVoxelsRefMut<'a, T, P>
where
    P: Into<usize>,
{
    pub fn get_pointed_voxel_info_mut(&mut self, ptr: P) -> &mut T {
        &mut self.infos[ptr.into()]
    }
}

impl<'a, T, P> GetExtent for ChunkVoxelsRefMut<'a, T, P> {
    fn get_extent(&self) -> &Extent {
        self.map.get_extent()
    }
}

impl<'a, T, P> GetLocalRefMut<T> for ChunkVoxelsRefMut<'a, T, P>
where
    P: Clone + Into<usize>,
{
    fn get_local_ref_mut(&mut self, p: &Point) -> &mut T {
        self.get_pointed_voxel_info_mut(self.map.get_local(p))
    }
}

/// A `ChunkedLatticeMap` that stores voxel data efficiently using palette compression.
#[derive(Deserialize, Serialize)]
pub struct ChunkedPaletteLatticeMap<T, P> {
    /// The palette of voxels that can be used in the map.
    pub infos: Vec<T>,
    /// Which voxels are used at specific points of the lattice.
    pub map: ChunkedLatticeMap<P>,
}

impl<T, P> ChunkedPaletteLatticeMap<T, P> {
    pub fn get_chunk(&self, chunk_key: &Point) -> Option<ChunkVoxelsRef<T, P>> {
        let ChunkedPaletteLatticeMap { map, infos } = self;

        map.get_chunk(chunk_key).map(|chunk_map| ChunkVoxelsRef {
            map: chunk_map,
            infos,
        })
    }

    pub fn get_chunk_mut(&mut self, chunk_key: &Point) -> Option<ChunkVoxelsRefMut<T, P>> {
        let ChunkedPaletteLatticeMap { map, infos } = self;

        map.get_mut_chunk(chunk_key)
            .map(move |chunk_map| ChunkVoxelsRefMut {
                map: chunk_map,
                infos,
            })
    }

    pub fn iter_chunks_ref(&self) -> impl Iterator<Item = (&Point, ChunkVoxelsRef<T, P>)> {
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

impl<T, P> ChunkedPaletteLatticeMap<T, P>
where
    P: Clone + Default,
{
    pub fn get_chunk_and_boundary(&self, chunk_key: &Point) -> LatticeVoxels<T, P> {
        let ChunkedPaletteLatticeMap { map, infos } = self;

        let chunk_and_boundary = map.get_chunk_and_boundary(chunk_key);

        LatticeVoxels {
            infos,
            map: chunk_and_boundary,
        }
    }

    pub fn iter_chunks_with_boundary(&self) -> impl Iterator<Item = (&Point, LatticeVoxels<T, P>)> {
        self.map
            .chunk_keys()
            .map(move |chunk_key| (chunk_key, self.get_chunk_and_boundary(chunk_key)))
    }
}

impl<T, P> MaybeGetWorldRef<T> for ChunkedPaletteLatticeMap<T, P>
where
    P: Clone + Into<usize>,
{
    fn maybe_get_world_ref(&self, p: &Point) -> Option<&T> {
        self.map
            .maybe_get_world(p)
            .map(|ptr| &self.infos[ptr.into()])
    }
}

impl<T, P> MaybeGetWorldRefMut<T> for ChunkedPaletteLatticeMap<T, P>
where
    P: Clone + Into<usize>,
{
    fn maybe_get_world_ref_mut(&mut self, p: &Point) -> Option<&mut T> {
        self.map
            .maybe_get_world(p)
            .map(move |ptr| &mut self.infos[ptr.into()])
    }
}
