use crate::{prelude::*, ChunkedLatticeMap, Extent, IsEmpty, Point, VecLatticeMap};
use serde::{Deserialize, Serialize};

/// One byte represents:
/// * a pointer to one of 127 possible voxel infos
/// * whether the voxel IsEmpty (for use in generic lattice algorithms)
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct VoxelInfoPtr {
    pub byte: u8,
}

const IS_EMPTY_MASK: u8 = 0b10000000;
const ADDRESS_MASK: u8 = !IS_EMPTY_MASK;

impl VoxelInfoPtr {
    pub fn new(address: u8, is_empty: bool) -> Self {
        let mut byte = address;
        if is_empty {
            byte |= IS_EMPTY_MASK;
        }

        VoxelInfoPtr { byte }
    }

    pub fn is_null(self) -> bool {
        self.address() == 0
    }

    pub fn address(self) -> usize {
        (self.byte & ADDRESS_MASK) as usize
    }
}

pub const NULL_VOXEL: VoxelInfoPtr = VoxelInfoPtr {
    byte: IS_EMPTY_MASK,
};

impl Default for VoxelInfoPtr {
    fn default() -> Self {
        NULL_VOXEL
    }
}

// Tell illatice3-mesh whether a voxel needs a cube mesh.
impl IsEmpty for VoxelInfoPtr {
    fn is_empty(&self) -> bool {
        (self.byte & IS_EMPTY_MASK) != 0
    }
}

/// An owned `VecLatticeMap` with a borrowed palette. The `VecLatticeMap` does not need to be a chunk.
pub struct LatticeVoxels<'a, T> {
    infos: &'a Vec<T>,
    pub lattice: VecLatticeMap<VoxelInfoPtr>,
}

impl<'a, T> LatticeVoxels<'a, T> {
    pub fn get_pointed_voxel_info(&'a self, ptr: VoxelInfoPtr) -> &'a T {
        &self.infos[ptr.address()]
    }
}

impl<'a, T> GetExtent for LatticeVoxels<'a, T> {
    fn get_extent(&self) -> Extent {
        self.lattice.get_extent()
    }
}

impl<'a, T> GetLocalRef<T> for LatticeVoxels<'a, T> {
    fn get_local_ref(&self, p: &Point) -> &T {
        self.get_pointed_voxel_info(self.lattice.get_local(p))
    }
}

/// A borrowed chunk and palette.
pub struct ChunkVoxelsRef<'a, T> {
    infos: &'a Vec<T>,
    pub lattice: &'a VecLatticeMap<VoxelInfoPtr>,
}

impl<'a, T> ChunkVoxelsRef<'a, T> {
    pub fn get_pointed_voxel_info(&self, ptr: VoxelInfoPtr) -> &T {
        &self.infos[ptr.address()]
    }
}

impl<'a, T> GetExtent for ChunkVoxelsRef<'a, T> {
    fn get_extent(&self) -> Extent {
        self.lattice.get_extent()
    }
}

impl<'a, T> GetLocalRef<T> for ChunkVoxelsRef<'a, T> {
    fn get_local_ref(&self, p: &Point) -> &T {
        self.get_pointed_voxel_info(self.lattice.get_local(p))
    }
}

/// A mutably borrowed chunk and palette.
pub struct ChunkVoxelsRefMut<'a, T> {
    infos: &'a mut Vec<T>,
    lattice: &'a mut VecLatticeMap<VoxelInfoPtr>,
}

impl<'a, T> ChunkVoxelsRefMut<'a, T> {
    pub fn get_pointed_voxel_info_mut(&mut self, ptr: VoxelInfoPtr) -> &mut T {
        &mut self.infos[ptr.address()]
    }
}

impl<'a, T> GetExtent for ChunkVoxelsRefMut<'a, T> {
    fn get_extent(&self) -> Extent {
        self.lattice.get_extent()
    }
}

impl<'a, T> GetLocalRefMut<T> for ChunkVoxelsRefMut<'a, T> {
    fn get_local_ref_mut(&mut self, p: &Point) -> &mut T {
        self.get_pointed_voxel_info_mut(self.lattice.get_local(p))
    }
}

/// A `ChunkedLatticeMap` that stores voxel data efficiently using palette compression.
#[derive(Deserialize, Serialize)]
pub struct ChunkedPaletteLatticeMap<T> {
    /// The palette of voxels that can be used in the lattice.
    pub infos: Vec<T>,
    /// Which voxels are used at specific points of the lattice.
    pub lattice: ChunkedLatticeMap<VoxelInfoPtr>,
}

impl<T> ChunkedPaletteLatticeMap<T> {
    pub fn get_chunk(&self, chunk_key: &Point) -> Option<ChunkVoxelsRef<T>> {
        let ChunkedPaletteLatticeMap { lattice, infos } = self;

        lattice
            .get_chunk(chunk_key)
            .map(|lattice| ChunkVoxelsRef { lattice, infos })
    }

    pub fn get_chunk_mut(&mut self, chunk_key: &Point) -> Option<ChunkVoxelsRefMut<T>> {
        let ChunkedPaletteLatticeMap { lattice, infos } = self;

        lattice
            .get_mut_chunk(chunk_key)
            .map(move |lattice| ChunkVoxelsRefMut { lattice, infos })
    }

    pub fn iter_chunks_ref(&self) -> impl Iterator<Item = (&Point, ChunkVoxelsRef<T>)> {
        self.lattice.iter_chunks().map(move |(chunk_key, lattice)| {
            (
                chunk_key,
                ChunkVoxelsRef {
                    lattice,
                    infos: &self.infos,
                },
            )
        })
    }

    pub fn get_chunk_and_boundary(&self, chunk_key: &Point) -> LatticeVoxels<T> {
        let ChunkedPaletteLatticeMap { lattice, infos } = self;

        let chunk_and_boundary = lattice.get_chunk_and_boundary(chunk_key);

        LatticeVoxels {
            infos,
            lattice: chunk_and_boundary,
        }
    }

    pub fn iter_chunks_with_boundary(&self) -> impl Iterator<Item = (&Point, LatticeVoxels<T>)> {
        self.lattice
            .chunk_keys()
            .map(move |chunk_key| (chunk_key, self.get_chunk_and_boundary(chunk_key)))
    }
}

impl<T> MaybeGetWorldRef<T> for ChunkedPaletteLatticeMap<T> {
    fn maybe_get_world_ref(&self, p: &Point) -> Option<&T> {
        self.lattice
            .maybe_get_world_ref(p)
            .map(|ptr| &self.infos[ptr.address()])
    }
}

impl<T> MaybeGetWorldRefMut<T> for ChunkedPaletteLatticeMap<T> {
    fn maybe_get_world_ref_mut(&mut self, p: &Point) -> Option<&mut T> {
        self.lattice
            .maybe_get_world_ref_mut(p)
            .cloned()
            .map(move |ptr| &mut self.infos[ptr.address()])
    }
}
