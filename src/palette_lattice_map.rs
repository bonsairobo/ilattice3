//! The various Get* traits are not implemented for these structures because it is customary for
//! users to store data in two places: the palette and the underlying lattice map. We don't
//! prescribe an implementation. Users may want to implement the Get* traits in multiple ways, and
//! they can use newtypes that wrap these types to do so.

use crate::{
    prelude::*, ChunkedLatticeMap, Extent, HasIndexer, Indexer, Point, VecLatticeMap,
    YLevelsIndexer,
};
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
    pub fn get_pointed_voxel_info(&'a self, ptr: P) -> &'a T {
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
    pub fn get_pointed_voxel_info(&self, ptr: P) -> &T {
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
    pub fn get_pointed_voxel_info_mut(&mut self, ptr: P) -> &mut T {
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
pub struct ChunkedPaletteLatticeMap<T, P, I = YLevelsIndexer> {
    /// The palette of voxels that can be used in the map.
    pub palette: Vec<T>,
    /// Which voxels are used at specific points of the lattice.
    pub map: ChunkedLatticeMap<P, I>,
}

impl<T, P, I> ChunkedPaletteLatticeMap<T, P, I> {
    pub fn get_chunk(&self, chunk_key: &Point) -> Option<ChunkVoxelsRef<T, P, I>> {
        let ChunkedPaletteLatticeMap { map, palette } = self;

        map.get_chunk(chunk_key).map(|chunk_map| ChunkVoxelsRef {
            map: chunk_map,
            palette,
        })
    }

    pub fn get_chunk_mut(&mut self, chunk_key: &Point) -> Option<ChunkVoxelsRefMut<T, P, I>> {
        let ChunkedPaletteLatticeMap { map, palette } = self;

        map.get_mut_chunk(chunk_key)
            .map(move |chunk_map| ChunkVoxelsRefMut {
                map: chunk_map,
                palette,
            })
    }

    pub fn iter_chunks_ref(&self) -> impl Iterator<Item = (&Point, ChunkVoxelsRef<T, P, I>)> {
        self.map.iter_chunks().map(move |(chunk_key, chunk_map)| {
            (
                chunk_key,
                ChunkVoxelsRef {
                    map: chunk_map,
                    palette: &self.palette,
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
