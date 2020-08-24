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
pub mod search;
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

pub use chunked_lattice_map::{Chunk, ChunkedLatticeMap};
pub use extent::{
    bounding_extent, copy_extent, copy_extent_to_position, fill_extent, map_extent, Extent,
};
pub use fn_lattice_map::FnLatticeMap;
pub use indexer::{PeriodicYLevelsIndexer, YLevelsIndexer};
pub use palette_lattice_map::{
    ChunkVoxelsRef, ChunkVoxelsRefMut, ChunkedPaletteLatticeMap, GetPaletteAddress, LatticeVoxels,
};
pub use point::Point;
pub use tile::Tile;
pub use transform::{Matrix, Transform};
pub use vec_lattice_map::VecLatticeMap;

pub trait IsEmpty {
    fn is_empty(&self) -> bool;
}

pub mod prelude {
    pub use crate::{
        extent::GetExtent,
        indexer::{HasIndexer, Indexer},
        lattice_map::{
            GetLinear, GetLinearRef, GetLinearRefMut, GetLocal, GetLocalRef, GetLocalRefMut,
            GetWorld, GetWorldBorrowable, GetWorldRef, GetWorldRefMut, LatticeMapIter,
            LatticeMapKeyValIterator, MaybeGetWorld, MaybeGetWorldRef, MaybeGetWorldRefMut,
        },
        Extent, IsEmpty, Point,
    };
}
