//! Data types, structures, and algorithms for 3D integer lattices (voxels)

// TODO: more robust linear algebra and transformations

mod chunked;
mod extent;
mod generic_lattice;
mod indexer;
mod lattice;
mod normal;
mod palette;
mod point;
mod tile;
mod transform;

#[cfg(feature = "vox")]
mod vox;

#[cfg(feature = "vox")]
pub use vox::{VoxColor, EMPTY_VOX_COLOR};

#[cfg(feature = "image")]
mod image;

#[cfg(test)]
mod test_util;

pub use chunked::{ChunkKeyIterator, ChunkedLattice, ChunkedLatticeIterator};
pub use extent::{bounding_extent, Extent, ExtentIterator};
pub use generic_lattice::{
    copy_extent, copy_extent_to_position, fill_extent, map_extent, GetExtent, GetLocal,
    GetLocalMut, GetWorld, GetWorldMut, MaybeGetWorld, MaybeGetWorldMut,
};
pub use indexer::{Indexer, PeriodicYLevelsIndexer, StatelessIndexer, YLevelsIndexer};
pub use lattice::{Lattice, LatticeKeyValIterator};
pub use normal::{
    closest_normal, normal_from_component_index, Direction, DirectionIndex, Normal, PlaneSpanInfo,
    ALL_DIRECTIONS, ALL_NORMALS,
};
pub use palette::{
    ChunkVoxelsRef, ChunkVoxelsRefMut, ChunkedPaletteLattice, LatticeVoxels, VoxelInfoPtr,
    NULL_VOXEL,
};
pub use point::Point;
pub use tile::Tile;
pub use transform::{Matrix, Transform, OCTAHEDRAL_GROUP, Z_STATIONARY_OCTAHEDRAL_GROUP};

pub trait IsEmpty: Eq + Sized {
    fn is_empty(&self) -> bool;
}

pub mod prelude {
    pub use crate::extent::{bounding_extent, Extent, ExtentIterator};
    pub use crate::generic_lattice::{
        copy_extent, copy_extent_to_position, fill_extent, map_extent, GetExtent, GetLocal,
        GetLocalMut, GetWorld, GetWorldMut, MaybeGetWorld, MaybeGetWorldMut,
    };
    pub use crate::lattice::{Lattice, LatticeKeyValIterator};
    pub use crate::point::Point;
}
