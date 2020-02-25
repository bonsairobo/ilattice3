//! Data types, structures, and algorithms for 3D integer lattices (voxels)

// TODO: more robust linear algebra and transformations

mod extent;
mod lattice;
mod normal;
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

pub use extent::{bounding_extent, Extent, ExtentIterator};
pub use lattice::{
    ChunkKeyIterator, ChunkedLattice, ChunkedLatticeIterator, Lattice, LatticeIndexer,
    LatticeKeyValIterator, PeriodicYLevelsIndexer, YLevelsIndexer,
};
pub use normal::{Direction, DirectionIndex, Normal, PlaneSpanInfo, ALL_DIRECTIONS, ALL_NORMALS};
pub use point::Point;
pub use tile::Tile;
pub use transform::{Transform, OCTAHEDRAL_GROUP, Z_STATIONARY_OCTAHEDRAL_GROUP};
