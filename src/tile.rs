use crate::{Extent, Indexer, Lattice, StatelessIndexer, YLevelsIndexer};

/// A container for voxels without any specified location, but they can be placed back into any
/// extent with the same size as their original extent, and their spatial order will be preserved.
/// This saves you from storing an `Extent` if you don't need it, and it also allows hashing and
/// comparison based only on the values.
#[derive(Clone, Eq, Hash, PartialEq)]
pub struct Tile<C, I = YLevelsIndexer> {
    data: Vec<C>,
    indexer: I,
}

impl<C, I: StatelessIndexer> Tile<C, I> {
    pub fn new(data: Vec<C>) -> Self {
        Tile {
            data,
            indexer: I::new(),
        }
    }
}

impl<C, I: Indexer> Tile<C, I> {
    pub fn new_with_indexer(indexer: I, data: Vec<C>) -> Self {
        Tile { data, indexer }
    }

    pub fn as_slice(&self) -> &[C] {
        &self.data
    }

    /// The primary constructor, copies the values in `extent` out of `lattice`, preserving the
    /// spatial structure.
    pub fn get_from_lattice<G: Clone + Into<C>>(
        lattice: &Lattice<G, I>,
        extent: &Extent,
    ) -> Tile<C, I> {
        Tile::new_with_indexer(
            lattice.get_indexer().clone(),
            lattice
                .serialize_extent(extent)
                .into_iter()
                .map(|g| g.into())
                .collect::<Vec<C>>(),
        )
    }

    /// Puts the tile in a specific location.
    pub fn put_in_extent(self, indexer: I, extent: Extent) -> Lattice<C, I> {
        Lattice::new_with_indexer(extent, indexer, self.data)
    }
}

impl<C: Clone, I: Clone + Indexer> Tile<C, I> {
    pub fn put_in_lattice<T: From<C>>(self, extent: &Extent, dst: &mut Lattice<T, I>) {
        let src = self.put_in_extent(dst.get_indexer().clone(), *extent);
        Lattice::copy_extent(&src, dst, extent);
    }
}
