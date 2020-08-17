//! It is customary for users to store data in two places: the palette and the underlying lattice
//! map. Users may want to implement the Get* traits in multiple ways, and they can use newtypes
//! that wrap these types to do so.

use crate::{prelude::*, ChunkedLatticeMap, Extent, Point, VecLatticeMap, YLevelsIndexer};
use serde::{Deserialize, Serialize};

pub trait GetPaletteAddress {
    fn get_palette_address(&self) -> usize;
}

/// An owned `VecLatticeMap` with a borrowed palette. The `VecLatticeMap` does not need to be a
/// chunk.
pub struct LatticeVoxels<'a, T, P, I = YLevelsIndexer> {
    pub palette: &'a Vec<T>,
    pub map: VecLatticeMap<P, I>,
}

impl<'a, T, P, I> LatticeVoxels<'a, T, P, I>
where
    P: GetPaletteAddress,
{
    pub fn get_pointed_voxel_info(&'a self, ptr: &P) -> &'a T {
        &self.palette[ptr.get_palette_address()]
    }
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

/// A borrowed chunk and palette.
pub struct ChunkVoxelsRef<'a, T, P, I> {
    pub palette: &'a Vec<T>,
    pub map: &'a VecLatticeMap<P, I>,
}

impl<'a, T, P, I> ChunkVoxelsRef<'a, T, P, I>
where
    P: GetPaletteAddress,
{
    pub fn get_pointed_voxel_info(&self, ptr: &P) -> &T {
        &self.palette[ptr.get_palette_address()]
    }
}

impl<'a, T, P, I> GetExtent for ChunkVoxelsRef<'a, T, P, I> {
    fn get_extent(&self) -> &Extent {
        self.map.get_extent()
    }
}

/// A mutably borrowed chunk and palette.
pub struct ChunkVoxelsRefMut<'a, T, P, I> {
    pub palette: &'a mut Vec<T>,
    pub map: &'a mut VecLatticeMap<P, I>,
}

impl<'a, T, P, I> ChunkVoxelsRefMut<'a, T, P, I>
where
    P: GetPaletteAddress,
{
    pub fn get_pointed_voxel_info_mut(&mut self, ptr: &P) -> &mut T {
        &mut self.palette[ptr.get_palette_address()]
    }
}

impl<'a, T, P, I> GetExtent for ChunkVoxelsRefMut<'a, T, P, I> {
    fn get_extent(&self) -> &Extent {
        self.map.get_extent()
    }
}

/// A `ChunkedLatticeMap` that stores voxel data efficiently using palette compression.
#[derive(Deserialize, Serialize)]
pub struct ChunkedPaletteLatticeMap<T, P, M = (), I = YLevelsIndexer> {
    /// The palette of voxels that can be used in the map.
    pub palette: Vec<T>,
    /// Which voxels are used at specific points of the lattice.
    pub map: ChunkedLatticeMap<P, M, I>,
}

impl<T, P, M, I> ChunkedPaletteLatticeMap<T, P, M, I> {
    pub fn get_chunk(&self, chunk_key: &Point) -> Option<ChunkVoxelsRef<T, P, I>> {
        let ChunkedPaletteLatticeMap { map, palette } = self;

        map.get_chunk(chunk_key).map(|chunk| ChunkVoxelsRef {
            map: &chunk.map,
            palette,
        })
    }

    pub fn get_chunk_mut(&mut self, chunk_key: &Point) -> Option<ChunkVoxelsRefMut<T, P, I>> {
        let ChunkedPaletteLatticeMap { map, palette } = self;

        map.get_mut_chunk(chunk_key)
            .map(move |chunk| ChunkVoxelsRefMut {
                map: &mut chunk.map,
                palette,
            })
    }

    pub fn iter_chunks_ref(&self) -> impl Iterator<Item = (&Point, ChunkVoxelsRef<T, P, I>)> {
        self.map.iter_chunks().map(move |(chunk_key, chunk)| {
            (
                chunk_key,
                ChunkVoxelsRef {
                    map: &chunk.map,
                    palette: &self.palette,
                },
            )
        })
    }
}

impl<T, P, M, I> ChunkedPaletteLatticeMap<T, P, M, I>
where
    P: Clone + Default,
    I: Indexer,
{
    pub fn copy_extent_into_new_map(&self, extent: Extent) -> LatticeVoxels<T, P, I> {
        let ChunkedPaletteLatticeMap { map, palette } = self;

        let voxels = map.copy_extent_into_new_map(extent);

        LatticeVoxels {
            palette,
            map: voxels,
        }
    }

    pub fn get_chunk_and_boundary(&self, chunk_key: &Point) -> LatticeVoxels<T, P, I> {
        let ChunkedPaletteLatticeMap { map, palette } = self;

        let chunk_and_boundary = map.get_chunk_and_boundary(chunk_key);

        LatticeVoxels {
            palette,
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
        self.get_pointed_voxel_info(self.map.get_linear_ref(i))
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
        let ptr = self.map.get_linear(i);

        self.get_pointed_voxel_info(&ptr)
    }
}

impl<'a, T, P, I> GetLocalRefMut for ChunkVoxelsRefMut<'a, T, P, I>
where
    P: Clone + GetPaletteAddress,
    I: Indexer,
{
    type Data = T;

    fn get_local_ref_mut(&mut self, p: &Point) -> &mut T {
        let ptr = self.map.get_local(p);

        self.get_pointed_voxel_info_mut(&ptr)
    }
}

impl<T, P, M, I> MaybeGetWorldRef for ChunkedPaletteLatticeMap<T, P, M, I>
where
    P: Clone + GetPaletteAddress,
    I: Indexer,
{
    type Data = T;

    fn maybe_get_world_ref(&self, p: &Point) -> Option<&T> {
        self.map
            .maybe_get_world(p)
            .map(|ptr| &self.palette[ptr.get_palette_address()])
    }
}

impl<T, P, M, I> MaybeGetWorldRefMut for ChunkedPaletteLatticeMap<T, P, M, I>
where
    P: Clone + GetPaletteAddress,
    I: Indexer,
{
    type Data = T;

    fn maybe_get_world_ref_mut(&mut self, p: &Point) -> Option<&mut T> {
        self.map
            .maybe_get_world(p)
            .map(move |ptr| &mut self.palette[ptr.get_palette_address()])
    }
}
