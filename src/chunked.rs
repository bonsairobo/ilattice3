use crate::{
    bounding_extent, lattice::LatticeKeyValIterator, Extent, Indexer, Lattice, Point, YLevelsIndexer
};

use std::collections::{hash_map, HashMap};

/// Stores a sparse lattice in chunks using a hash map.
pub struct ChunkedLattice<T> {
    chunk_size: Point,
    map: HashMap<Point, Lattice<T>>,
}

impl<T> ChunkedLattice<T> {
    pub fn new(chunk_size: Point) -> Self {
        assert!(chunk_size > [0, 0, 0].into());

        Self {
            chunk_size,
            map: HashMap::new(),
        }
    }

    pub fn get_extent(&self) -> Extent {
        assert!(!self.map.is_empty());
        let extrema_iter = self
            .map
            .values()
            .map(|lat| lat.get_extent().get_world_max())
            .chain(self.map.values().map(|lat| lat.get_extent().get_minimum()));

        bounding_extent(extrema_iter)
    }

    pub fn chunk_key(&self, point: &Point) -> Point {
        *point / self.chunk_size
    }

    fn key_extent(&self, extent: &Extent) -> Extent {
        let key_min = self.chunk_key(&extent.get_minimum());
        let key_max = self.chunk_key(&extent.get_world_max());

        Extent::from_min_and_world_max(key_min, key_max)
    }

    pub fn extent_for_chunk_key(&self, key: &Point) -> Extent {
        let min = *key * self.chunk_size;
        let local_sup = self.chunk_size;

        Extent::from_min_and_local_supremum(min, local_sup)
    }

    pub fn get_chunk(&self, key: &Point) -> Option<&Lattice<T>> {
        self.map.get(key)
    }

    pub fn get_mut_chunk(&mut self, key: &Point) -> Option<&mut Lattice<T>> {
        self.map.get_mut(key)
    }

    pub fn get_chunk_containing_point(&self, point: &Point) -> Option<&Lattice<T>> {
        self.get_chunk(&self.chunk_key(point))
    }

    pub fn get_mut_chunk_containing_point(&mut self, point: &Point) -> Option<&mut Lattice<T>> {
        self.get_mut_chunk(&self.chunk_key(point))
    }

    /// Returns an iterator over all points and corresponding values in the given extent. If chunks
    /// are missing from the extent, then their points will not be yielded.
    pub fn iter_point_values(
        &self,
        extent: Extent,
    ) -> ChunkedLatticeIterator<T, std::vec::IntoIter<&Lattice<T, YLevelsIndexer>>, YLevelsIndexer>
    {
        let lattices = self
            .key_extent(&extent)
            .into_iter()
            .filter_map(|k| self.map.get(&k))
            .collect::<Vec<&Lattice<T>>>()
            .into_iter();

        ChunkedLatticeIterator {
            full_extent: extent,
            lattices,
            lattice_iter: None,
        }
    }

    /// An iterator over the chunk keys (sparse points in the space of chunk coordinates).
    pub fn chunk_keys(&self) -> ChunkKeyIterator<T> {
        ChunkKeyIterator {
            map_key_iter: self.map.keys(),
        }
    }

    /// Get the data for point `p` if `p` exists.
    pub fn get_world(&self, p: &Point) -> Option<&T> {
        self.get_chunk_containing_point(p)
            .map(|chunk| chunk.get_world(p))
    }

    /// Get mutable data for point `p` if `p` exists.
    pub fn get_mut_world(&mut self, p: &Point) -> Option<&mut T> {
        self.get_mut_chunk_containing_point(p)
            .map(|chunk| chunk.get_mut_world(p))
    }
}

/// An iterator over the chunk keys (sparse points in the space of chunk coordinates).
pub struct ChunkKeyIterator<'a, T> {
    map_key_iter: hash_map::Keys<'a, Point, Lattice<T>>,
}

impl<'a, T> Iterator for ChunkKeyIterator<'a, T> {
    type Item = &'a Point;

    fn next(&mut self) -> Option<Self::Item> {
        self.map_key_iter.next()
    }
}

impl<T: Clone + Default> ChunkedLattice<T> {
    pub fn copy_lattice_into_chunks(&mut self, lattice: &Lattice<T>, fill_default: T) {
        for key in self.key_extent(&lattice.get_extent()) {
            let chunk_extent = self.extent_for_chunk_key(&key);
            let chunk = self
                .map
                .entry(key)
                .or_insert_with(|| Lattice::fill(chunk_extent, fill_default.clone()));
            Lattice::copy_extent(
                lattice,
                chunk,
                &chunk_extent.intersection(&lattice.get_extent()),
            );
        }
    }

    /// `fill_default` will be the value of points outside the given extent that nonetheless must be
    /// filled in each non-sparse chunk.
    pub fn fill_extent(&mut self, extent: &Extent, val: T, fill_default: T) {
        let fill_lat = Lattice::fill(*extent, val);
        self.copy_lattice_into_chunks(&fill_lat, fill_default);
    }

    pub fn copy_extent_into_new_lattice(&self, extent: Extent) -> Lattice<T> {
        let mut new_lattice = Lattice::fill(extent, T::default());
        for (p, val) in self.iter_point_values(extent) {
            *new_lattice.get_mut_world(&p) = val;
        }

        new_lattice
    }

    pub fn copy_into_new_lattice(&self) -> Lattice<T> {
        self.copy_extent_into_new_lattice(self.get_extent())
    }
}

/// An iterator over the points in a `ChunkedLattice`.
pub struct ChunkedLatticeIterator<'a, T, It, In>
where
    It: Iterator<Item = &'a Lattice<T, In>>,
{
    full_extent: Extent,
    lattices: It,
    lattice_iter: Option<LatticeKeyValIterator<'a, T, In>>,
}

impl<'a, T, It, In> ChunkedLatticeIterator<'a, T, It, In>
where
    T: Clone,
    It: Iterator<Item = &'a Lattice<T, In>>,
    In: Indexer,
{
    fn move_to_next_lattice(&mut self) {
        self.lattice_iter = self.lattices.next().map(|l| {
            let extent_iter = l.get_extent().intersection(&self.full_extent).into_iter();

            LatticeKeyValIterator::new(l, extent_iter)
        });
    }
}

impl<'a, T, It, In> Iterator for ChunkedLatticeIterator<'a, T, It, In>
where
    T: Clone,
    It: Iterator<Item = &'a Lattice<T, In>>,
    In: Indexer,
{
    type Item = (Point, T);

    fn next(&mut self) -> Option<Self::Item> {
        if self.lattice_iter.is_none() {
            self.move_to_next_lattice();
        }

        let mut next = None;
        while let Some(i) = &mut self.lattice_iter {
            let n = i.next();
            if n.is_none() {
                self.move_to_next_lattice();
            } else {
                next = n;
                break;
            }
        }

        next
    }
}

// ████████╗███████╗███████╗████████╗███████╗
// ╚══██╔══╝██╔════╝██╔════╝╚══██╔══╝██╔════╝
//    ██║   █████╗  ███████╗   ██║   ███████╗
//    ██║   ██╔══╝  ╚════██║   ██║   ╚════██║
//    ██║   ███████╗███████║   ██║   ███████║
//    ╚═╝   ╚══════╝╚══════╝   ╚═╝   ╚══════╝

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::assert_elements_eq;

    #[test]
    fn test_chunked_lattice_iterator() {
        // All points less than [0, 0, 0] in the query extent will have value 1.
        let sublattice1 = Lattice::fill(Extent::from_center_and_radius([-7, -7, -7].into(), 7), 1);
        // Only [1, 1, 1] with have value 2 in the query extent.
        let sublattice2 = Lattice::fill(Extent::from_center_and_radius([8, 8, 8].into(), 7), 2);
        // All other point values will not change (stay 0).

        let mut chunked_lattice: ChunkedLattice<u32> = ChunkedLattice::new([4, 4, 4].into());
        chunked_lattice.copy_lattice_into_chunks(&sublattice1, 0);
        chunked_lattice.copy_lattice_into_chunks(&sublattice2, 0);

        let query_extent = Extent::from_center_and_radius([0, 0, 0].into(), 1);
        let points: Vec<_> = chunked_lattice.iter_point_values(query_extent).collect();

        assert_elements_eq(
            &points,
            &vec![
                ([-1, -1, -1].into(), 1),
                ([-1, -1, 0].into(), 1),
                ([-1, -1, 1].into(), 0),
                ([-1, 0, -1].into(), 1),
                ([-1, 0, 0].into(), 1),
                ([-1, 0, 1].into(), 0),
                ([-1, 1, -1].into(), 0),
                ([-1, 1, 0].into(), 0),
                ([-1, 1, 1].into(), 0),
                ([0, -1, -1].into(), 1),
                ([0, -1, 0].into(), 1),
                ([0, -1, 1].into(), 0),
                ([0, 0, -1].into(), 1),
                ([0, 0, 0].into(), 1),
                ([0, 0, 1].into(), 0),
                ([0, 1, -1].into(), 0),
                ([0, 1, 0].into(), 0),
                ([0, 1, 1].into(), 0),
                ([1, -1, -1].into(), 0),
                ([1, -1, 0].into(), 0),
                ([1, -1, 1].into(), 0),
                ([1, 0, -1].into(), 0),
                ([1, 0, 0].into(), 0),
                ([1, 0, 1].into(), 0),
                ([1, 1, -1].into(), 0),
                ([1, 1, 0].into(), 0),
                ([1, 1, 1].into(), 2),
            ],
        );
    }
}
