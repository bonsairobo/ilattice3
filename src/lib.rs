//! Data types, structures, and algorithms for 3D integer lattices (voxels)

mod extent;
mod lattice;
mod normal;
mod point;

#[cfg(test)]
mod test_util;

pub use extent::{bounding_extent, Extent, ExtentIterator};
pub use lattice::{
    ChunkKeyIterator, ChunkedLattice, ChunkedLatticeIterator, Lattice, LatticeIndexer,
    LatticeKeyValIterator, PeriodicYLevelsIndexer,
};
pub use normal::{Direction, DirectionIndex, Normal, PlaneSpanInfo, ALL_DIRECTIONS, ALL_NORMALS};
pub use point::Point;
