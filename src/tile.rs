use crate::{copy_extent, Extent, Indexer, Point, VecLatticeMap, YLevelsIndexer};

use std::marker::PhantomData;

/// A container for voxels without any specified location, but they can be placed back into any
/// extent with the same size as their original extent, and their spatial order will be preserved.
/// This allows hashing and comparison based only on the values and dimensions.
#[derive(Clone, Eq, Hash, PartialEq)]
pub struct Tile<C, I = YLevelsIndexer> {
    data: Vec<C>,
    dimensions: Point,
    indexer: PhantomData<I>,
}

impl<C, I: Indexer> Tile<C, I> {
    pub fn new(data: Vec<C>, dimensions: Point) -> Self {
        Tile {
            data,
            dimensions,
            indexer: Default::default(),
        }
    }
}

impl<C, I: Indexer> Tile<C, I> {
    pub fn as_slice(&self) -> &[C] {
        &self.data
    }

    /// The primary constructor, copies the values in `extent` out of `lattice`, preserving the
    /// spatial structure.
    pub fn get_from_map<G: Clone + Into<C>>(
        map: &VecLatticeMap<G, I>,
        extent: &Extent,
    ) -> Tile<C, I> {
        Tile::new(
            map.serialize_extent(extent)
                .into_iter()
                .map(|g| g.into())
                .collect::<Vec<C>>(),
            *extent.get_local_supremum(),
        )
    }

    /// Puts the tile in a specific location.
    pub fn put_in_extent(self, extent: Extent) -> VecLatticeMap<C, I> {
        assert_eq!(extent.get_local_supremum(), &self.dimensions);

        VecLatticeMap::new(extent, self.data)
    }
}

impl<C: Clone, I: Clone + Indexer> Tile<C, I> {
    pub fn put_in_map<T: From<C>>(self, extent: &Extent, dst: &mut VecLatticeMap<T, I>) {
        let src = self.put_in_extent(*extent);
        copy_extent(&src, dst, extent);
    }
}
