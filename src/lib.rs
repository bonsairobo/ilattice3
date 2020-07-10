//! Data types, structures, and algorithms for 3D integer lattices (voxels)

pub mod chunked_lattice_map;
pub mod extent;
pub mod fn_lattice_map;
pub mod indexer;
pub mod lattice_map;
pub mod normal;
pub mod palette_lattice_map;
pub mod point;
pub mod tile;
pub mod transform;
pub mod vec_lattice_map;

#[cfg(feature = "vox")]
mod vox;

#[cfg(feature = "vox")]
pub use vox::{VoxColor, EMPTY_VOX_COLOR};

#[cfg(feature = "image")]
mod image;

#[cfg(test)]
mod test_util;

pub use chunked_lattice_map::{ChunkKeyIterator, ChunkedLatticeMap, ChunkedLatticeMapIterator};
pub use extent::{
    bounding_extent, copy_extent, copy_extent_to_position, fill_extent, map_extent, Extent,
    ExtentIterator, GetExtent,
};
pub use fn_lattice_map::FnLatticeMap;
pub use indexer::{Indexer, PeriodicYLevelsIndexer, StatelessIndexer, YLevelsIndexer};
pub use lattice_map::{
    GetLocal, GetLocalRef, GetLocalRefMut, GetWorld, GetWorldBorrowable, GetWorldRef,
    GetWorldRefMut, LatticeMapKeyValIterator, MaybeGetWorld, MaybeGetWorldRef, MaybeGetWorldRefMut,
};
pub use normal::{
    closest_normal, normal_from_component_index, Direction, DirectionIndex, Normal, PlaneSpanInfo,
    ALL_DIRECTIONS, ALL_NORMALS,
};
pub use palette_lattice_map::{
    ChunkVoxelsRef, ChunkVoxelsRefMut, ChunkedPaletteLatticeMap, LatticeVoxels, VoxelInfoPtr,
    NULL_VOXEL,
};
pub use point::Point;
pub use tile::Tile;
pub use transform::{Matrix, Transform, OCTAHEDRAL_GROUP, Z_STATIONARY_OCTAHEDRAL_GROUP};
pub use vec_lattice_map::VecLatticeMap;

pub trait IsEmpty: Eq + Sized {
    fn is_empty(&self) -> bool;
}

pub mod prelude {
    pub use crate::extent::GetExtent;
    pub use crate::lattice_map::{
        GetLocal, GetLocalRef, GetLocalRefMut, GetWorld, GetWorldBorrowable, GetWorldRef,
        GetWorldRefMut, MaybeGetWorld, MaybeGetWorldRef, MaybeGetWorldRefMut,
    };
    pub use crate::point::Point;
}
