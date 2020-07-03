use crate::{ChunkedLattice, IsEmpty, Lattice, Point};
use serde::{Deserialize, Serialize};

/// One byte represents:
/// * a pointer to one of 127 possible voxel infos
/// * whether the voxel IsEmpty (for use in generic lattice algorithms)
#[derive(
    Clone, Copy, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize,
)]
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

// Tell illatice3-mesh whether a voxel needs a cube mesh.
impl IsEmpty for VoxelInfoPtr {
    fn is_empty(&self) -> bool {
        (self.byte & IS_EMPTY_MASK) != 0
    }
}

pub struct ChunkVoxelsMut<'a, T> {
    /// The palette of voxels that can be used in the lattice.
    infos: &'a mut Vec<T>,
    /// Which voxels are used at specific points of the lattice.
    lattice: &'a mut Lattice<VoxelInfoPtr>,
}

impl<'a, T> ChunkVoxelsMut<'a, T> {
    pub fn get_voxel_info_mut(&'a mut self, point: &Point) -> &'a mut T {
        &mut self.infos[self.lattice.get_world(point).address()]
    }
}

pub struct ChunkVoxels<'a, T> {
    /// The palette of voxels that can be used in the lattice.
    infos: &'a Vec<T>,
    /// Which voxels are used at specific points of the lattice.
    pub lattice: &'a Lattice<VoxelInfoPtr>,
}

impl<'a, T> ChunkVoxels<'a, T> {
    pub fn get_pointed_voxel_info(&'a self, ptr: VoxelInfoPtr) -> &'a T {
        &self.infos[ptr.address()]
    }

    pub fn get_voxel_info(&'a self, point: &Point) -> &'a T {
        self.get_pointed_voxel_info(*self.lattice.get_world(point))
    }
}

/// A `ChunkedLattice` that stores voxel data efficiently using palette compression.
#[derive(Deserialize, Serialize)]
pub struct ChunkedPaletteLattice<T> {
    /// The palette of voxels that can be used in the lattice.
    pub infos: Vec<T>,
    /// Which voxels are used at specific points of the lattice.
    pub lattice: ChunkedLattice<VoxelInfoPtr>,
}

impl<T> ChunkedPaletteLattice<T> {
    pub fn get_voxel_info(&self, point: &Point) -> Option<&T> {
        self.lattice
            .get_world(point)
            .cloned()
            .map(move |ptr| &self.infos[ptr.address()])
    }

    pub fn get_voxel_info_mut(&mut self, point: &Point) -> Option<&mut T> {
        self.lattice
            .get_world(point)
            .cloned()
            .map(move |ptr| &mut self.infos[ptr.address()])
    }

    pub fn get_chunk(&self, chunk_key: &Point) -> Option<ChunkVoxels<T>> {
        let ChunkedPaletteLattice { lattice, infos } = self;

        lattice
            .get_chunk(chunk_key)
            .map(|lattice| ChunkVoxels { lattice, infos })
    }

    pub fn get_chunk_mut(&mut self, chunk_key: &Point) -> Option<ChunkVoxelsMut<T>> {
        let ChunkedPaletteLattice { lattice, infos } = self;

        lattice
            .get_mut_chunk(chunk_key)
            .map(move |lattice| ChunkVoxelsMut { lattice, infos })
    }

    pub fn iter_chunks(&self) -> impl Iterator<Item = (&Point, ChunkVoxels<T>)> {
        self.lattice.iter_chunks().map(move |(chunk_key, lattice)| {
            (
                chunk_key,
                ChunkVoxels {
                    lattice,
                    infos: &self.infos,
                },
            )
        })
    }
}
