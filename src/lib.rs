//! Data types, structures, and algorithms for 3D integer lattices (voxels)

pub mod algos;
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

pub use algos::{find_surface_voxels, flood_fill_local, flood_fill_world};
pub use chunked_lattice_map::{ChunkKeyIterator, ChunkedLatticeMap, ChunkedLatticeMapIterator};
pub use extent::{
    bounding_extent, copy_extent, copy_extent_to_position, fill_extent, map_extent, Extent,
    ExtentIterator, GetExtent,
};
pub use fn_lattice_map::FnLatticeMap;
pub use indexer::{HasIndexer, Indexer, PeriodicYLevelsIndexer, YLevelsIndexer};
pub use lattice_map::{
    GetLinear, GetLinearRef, GetLinearRefMut, GetLocal, GetLocalRef, GetLocalRefMut, GetWorld,
    GetWorldBorrowable, GetWorldRef, GetWorldRefMut, LatticeMapIter, LatticeMapKeyValIterator,
    MaybeGetWorld, MaybeGetWorldRef, MaybeGetWorldRefMut,
};
pub use normal::{
    closest_normal, normal_from_component_index, Direction, DirectionIndex, Normal, PlaneSpanInfo,
    ALL_DIRECTIONS, ALL_NORMALS,
};
pub use palette_lattice_map::{
    ChunkVoxelsRef, ChunkVoxelsRefMut, ChunkedPaletteLatticeMap, GetPaletteAddress, LatticeVoxels,
};
pub use point::Point;
pub use tile::Tile;
pub use transform::{Matrix, Transform, OCTAHEDRAL_GROUP, Z_STATIONARY_OCTAHEDRAL_GROUP};
pub use vec_lattice_map::VecLatticeMap;

pub trait IsEmpty {
    fn is_empty(&self) -> bool;
}

pub mod prelude {
    pub use crate::extent::GetExtent;
    pub use crate::indexer::HasIndexer;
    pub use crate::lattice_map::{
        GetLinear, GetLinearRef, GetLinearRefMut, GetLocal, GetLocalRef, GetLocalRefMut, GetWorld,
        GetWorldBorrowable, GetWorldRef, GetWorldRefMut, LatticeMapIter, LatticeMapKeyValIterator,
        MaybeGetWorld, MaybeGetWorldRef, MaybeGetWorldRefMut,
    };
    pub use crate::point::Point;
}

pub const FACE_ADJACENT_OFFSETS: [Point; 6] = [
    Point { x: 1, y: 0, z: 0 },
    Point { x: 0, y: 1, z: 0 },
    Point { x: 0, y: 0, z: 1 },
    Point { x: -1, y: 0, z: 0 },
    Point { x: 0, y: -1, z: 0 },
    Point { x: 0, y: 0, z: -1 },
];

#[rustfmt::skip]
pub const ALL_ADJACENT_OFFSETS: [Point; 26] = [
    Point { x: -1, y: -1, z: -1 },
    Point { x: -1, y: -1, z:  0 },
    Point { x: -1, y: -1, z:  1 },
    Point { x: -1, y:  0, z: -1 },
    Point { x: -1, y:  0, z:  0 },
    Point { x: -1, y:  0, z:  1 },
    Point { x: -1, y:  1, z: -1 },
    Point { x: -1, y:  1, z:  0 },
    Point { x: -1, y:  1, z:  1 },
    Point { x:  0, y: -1, z: -1 },
    Point { x:  0, y: -1, z:  0 },
    Point { x:  0, y: -1, z:  1 },
    Point { x:  0, y:  0, z: -1 },
    Point { x:  0, y:  0, z:  1 },
    Point { x:  0, y:  1, z: -1 },
    Point { x:  0, y:  1, z:  0 },
    Point { x:  0, y:  1, z:  1 },
    Point { x:  1, y: -1, z: -1 },
    Point { x:  1, y: -1, z:  0 },
    Point { x:  1, y: -1, z:  1 },
    Point { x:  1, y:  0, z: -1 },
    Point { x:  1, y:  0, z:  0 },
    Point { x:  1, y:  0, z:  1 },
    Point { x:  1, y:  1, z: -1 },
    Point { x:  1, y:  1, z:  0 },
    Point { x:  1, y:  1, z:  1 },
];

pub const CUBE_CORNERS: [Point; 8] = [
    Point { x: 0, y: 0, z: 0 },
    Point { x: 1, y: 0, z: 0 },
    Point { x: 0, y: 1, z: 0 },
    Point { x: 1, y: 1, z: 0 },
    Point { x: 0, y: 0, z: 1 },
    Point { x: 1, y: 0, z: 1 },
    Point { x: 0, y: 1, z: 1 },
    Point { x: 1, y: 1, z: 1 },
];
