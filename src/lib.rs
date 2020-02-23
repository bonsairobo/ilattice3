//! Data types, structures, and algorithms for 3D integer lattices (voxels)

mod extent;
mod lattice;
mod normal;
mod point;

#[cfg(feature = "vox")]
mod vox;

#[cfg(feature = "vox")]
pub use vox::{VoxColor, EMPTY_VOX_COLOR};

#[cfg(test)]
mod test_util;

pub use extent::{bounding_extent, Extent, ExtentIterator};
pub use lattice::{
    ChunkKeyIterator, ChunkedLattice, ChunkedLatticeIterator, Lattice, LatticeIndexer,
    LatticeKeyValIterator, PeriodicYLevelsIndexer, YLevelsIndexer,
};
pub use normal::{Direction, DirectionIndex, Normal, PlaneSpanInfo, ALL_DIRECTIONS, ALL_NORMALS};
pub use point::Point;
