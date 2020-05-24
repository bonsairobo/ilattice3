//! Data types, structures, and algorithms for 3D integer lattices (voxels)

// TODO: more robust linear algebra and transformations

mod chunked;
mod extent;
mod indexer;
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

pub use chunked::{ChunkKeyIterator, ChunkedLattice, ChunkedLatticeIterator};
pub use extent::{bounding_extent, Extent, ExtentIterator};
pub use indexer::{Indexer, PeriodicYLevelsIndexer, StatelessIndexer, YLevelsIndexer};
pub use lattice::{Lattice, LatticeKeyValIterator};
pub use normal::{
    normal_from_component_index, Direction, DirectionIndex, Normal, PlaneSpanInfo, ALL_DIRECTIONS,
    ALL_NORMALS
};
pub use point::Point;
pub use tile::Tile;
pub use transform::{Matrix, Transform, OCTAHEDRAL_GROUP, Z_STATIONARY_OCTAHEDRAL_GROUP};

pub trait IsEmpty: Eq + Sized {
    fn is_empty(&self) -> bool;
}
