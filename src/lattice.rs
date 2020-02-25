use crate::{
    extent::{bounding_extent, Extent, ExtentIterator},
    transform::Transform,
    Point,
};

use std::collections::{hash_map, HashMap};

pub trait PeriodicLatticeIndexer: LatticeIndexer {}

pub trait LatticeIndexer {
    /// `s` is the local strict supremum of an extent. `p` is a local point.
    fn index_from_local_point(s: &Point, p: &Point) -> usize;

    /// `s` is the local strict supremum of an extent. `index` is a linear index.
    fn local_point_from_index(s: &Point, index: usize) -> Point;
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct YLevelsIndexer;

impl LatticeIndexer for YLevelsIndexer {
    fn index_from_local_point(s: &Point, p: &Point) -> usize {
        // This scheme is chosen for ease of hand-crafting voxel maps.
        (p.y * s.x * s.z + p.z * s.x + p.x) as usize
    }

    fn local_point_from_index(s: &Point, index: usize) -> Point {
        assert!(index <= std::i32::MAX as usize);
        let index = index as i32;
        let xz_area = s.x * s.z;
        let y = index / xz_area;
        let rem = index - y * xz_area;
        let z = rem / s.x;
        let rem = rem - z * s.x;
        let x = rem;

        Point::new(x, y, z)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PeriodicYLevelsIndexer;

impl PeriodicLatticeIndexer for PeriodicYLevelsIndexer {}

impl LatticeIndexer for PeriodicYLevelsIndexer {
    fn index_from_local_point(s: &Point, p: &Point) -> usize {
        let px = p.x.rem_euclid(s.x);
        let py = p.y.rem_euclid(s.y);
        let pz = p.z.rem_euclid(s.z);

        (py * s.x * s.z + pz * s.x + px) as usize
    }

    /// Returns the canonical point which is contained in the extent.
    fn local_point_from_index(s: &Point, index: usize) -> Point {
        YLevelsIndexer::local_point_from_index(s, index)
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Lattice<T, I = YLevelsIndexer> {
    pub indexer: I,
    extent: Extent,
    values: Vec<T>,
}

impl<T: Clone, I: Clone + LatticeIndexer> Lattice<T, I> {
    /// Map every point by `tfm`. This function will assert `tfm.is_octahedral` in debug mode.
    pub fn apply_octahedral_transform(&self, tfm: &Transform) -> Self {
        debug_assert!(tfm.is_octahedral());

        let extent = self.get_extent();
        let volume = extent.volume();

        let tfm_extent = tfm.apply_to_extent(&extent);

        let mut new_values = Vec::with_capacity(volume);
        unsafe {
            new_values.set_len(volume);
        }
        let mut tfm_lattice = Self::new_with_indexer(tfm_extent, self.indexer.clone(), new_values);

        // PERF: this is not the most efficient, but it is very simple.
        for p in extent {
            let tfm_p = tfm.apply_to_point(&p);
            *tfm_lattice.get_mut_world(&tfm_p) = self.get_world(&p).clone();
        }

        tfm_lattice
    }
}

impl<T: Clone, I: LatticeIndexer> Lattice<T, I> {
    pub fn fill_with_indexer(indexer: I, extent: Extent, init_val: T) -> Self {
        Lattice {
            extent,
            indexer,
            values: vec![init_val; extent.volume()],
        }
    }

    /// Returns a vec of the data in `extent`, ordered linearly by `I: LatticeIndexer`. A lattice
    /// can be recreated from the vec using `Lattice::<T, I>::new_with_indexer`.
    pub fn serialize_extent(&self, extent: &Extent) -> Vec<T> {
        let num_elements = extent.volume();
        let mut data = Vec::with_capacity(num_elements);
        unsafe {
            data.set_len(num_elements);
        }
        for p in extent {
            let i = I::index_from_local_point(
                extent.get_local_supremum(),
                &extent.local_point_from_world_point(&p),
            );
            data[i] = self.get_world(&p).clone();
        }

        data
    }

    pub fn fill_extent(&mut self, extent: &Extent, val: T) {
        for p in extent {
            *self.get_mut_world(&p) = val.clone();
        }
    }

    pub fn copy_extent_to_position<S: From<T>, J: LatticeIndexer>(
        src: &Self,
        dst: &mut Lattice<S, J>,
        dst_position: &Point,
        extent: &Extent,
    ) {
        for p in extent {
            let p_dst = *dst_position + p - extent.get_minimum();
            *dst.get_mut_world(&p_dst) = src.get_world(&p).clone().into();
        }
    }

    pub fn copy_extent<S: From<T>, J: LatticeIndexer>(
        src: &Self, dst: &mut Lattice<S, J>, extent: &Extent
    ) {
        Self::copy_extent_to_position(src, dst, &extent.get_minimum(), extent)
    }

    pub fn map_extent<S, F, J>(src: &Self, dst: &mut Lattice<S, J>, extent: &Extent, f: F)
    where
        F: Fn(&T) -> S,
        J: LatticeIndexer,
    {
        for p in extent {
            *dst.get_mut_world(&p) = f(src.get_world(&p));
        }
    }
}

impl<T: Clone> Lattice<T> {
    pub fn fill(extent: Extent, init_val: T) -> Self {
        Self::fill_with_indexer(YLevelsIndexer {}, extent, init_val)
    }
}

impl<T: Clone + Default> Lattice<T> {
    pub fn copy_extent_into_new_lattice(&self, extent: &Extent) -> Self {
        let mut copy = Lattice::fill(*extent, T::default());
        Self::copy_extent(self, &mut copy, extent);

        copy
    }
}

impl<T> Lattice<T> {
    pub fn new(extent: Extent, values: Vec<T>) -> Self {
        Lattice {
            extent,
            indexer: YLevelsIndexer {},
            values,
        }
    }

    pub fn new_at_origin(sup: Point, values: Vec<T>) -> Self {
        let extent = Extent::from_min_and_world_supremum([0, 0, 0].into(), sup);

        Self::new(extent, values)
    }
}

impl<T, I: Clone + LatticeIndexer> Lattice<T, I> {
    pub fn map<F, S>(&self, f: F) -> Lattice<S, I>
    where
        F: Fn(&T) -> S,
    {
        Lattice::new_with_indexer(
            self.get_extent(),
            self.indexer.clone(),
            self.values.iter().map(f).collect(),
        )
    }
}

impl<T, I: LatticeIndexer> Lattice<T, I> {
    pub fn new_with_indexer(extent: Extent, indexer: I, values: Vec<T>) -> Self {
        Lattice {
            extent,
            indexer,
            values,
        }
    }

    pub fn get_extent(&self) -> Extent {
        self.extent
    }

    pub fn get_linear(&self, index: usize) -> &T {
        &self.values[index]
    }

    pub fn get_mut_linear(&mut self, index: usize) -> &mut T {
        &mut self.values[index]
    }

    pub fn index_from_local_point(&self, p: &Point) -> usize {
        let local_sup = self.extent.get_local_supremum();

        I::index_from_local_point(&local_sup, p)
    }

    pub fn local_point_from_index(&self, index: usize) -> Point {
        let local_sup = self.extent.get_local_supremum();

        I::local_point_from_index(&local_sup, index)
    }

    pub fn get_local(&self, p: &Point) -> &T {
        self.get_linear(self.index_from_local_point(p))
    }

    pub fn get_mut_local(&mut self, p: &Point) -> &mut T {
        let i = self.index_from_local_point(p);

        &mut self.values[i]
    }

    pub fn index_from_world_point(&self, p: &Point) -> usize {
        self.index_from_local_point(&self.extent.local_point_from_world_point(p))
    }

    pub fn get_world(&self, p: &Point) -> &T {
        self.get_local(&self.extent.local_point_from_world_point(p))
    }

    pub fn get_mut_world(&mut self, p: &Point) -> &mut T {
        self.get_mut_local(&self.extent.local_point_from_world_point(p))
    }

    pub fn all<F>(&self, extent: &Extent, f: F) -> bool
    where
        F: Fn(&T) -> bool,
    {
        for p in extent {
            if !f(self.get_world(&p)) {
                return false;
            }
        }

        true
    }

    pub fn some<F>(&self, extent: &Extent, f: F) -> bool
    where
        F: Fn(&T) -> bool,
    {
        for p in extent {
            if f(self.get_world(&p)) {
                return true;
            }
        }

        false
    }
}

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
            .map(|lat| lat.extent.get_world_max())
            .chain(self.map.values().map(|lat| lat.extent.get_minimum()));

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

    pub fn get_chunk_containing_point(&self, point: &Point) -> Option<&Lattice<T>> {
        self.map.get(&self.chunk_key(point))
    }

    /// Returns an iterator over all points and corresponding values in the given extent. If chunks
    /// are missing from the extent, then their points will not be yielded.
    pub fn iter_point_values(
        &self,
        extent: Extent,
    ) -> ChunkedLatticeIterator<T, std::vec::IntoIter<&Lattice<T>>> {
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

    pub fn chunk_keys(&self) -> ChunkKeyIterator<T> {
        ChunkKeyIterator {
            map_key_iter: self.map.keys(),
        }
    }

    pub fn get_world(&self, p: &Point) -> Option<&T> {
        self.get_chunk_containing_point(p)
            .map(|chunk| chunk.get_world(p))
    }
}

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

pub struct LatticeKeyValIterator<'a, T> {
    lattice: &'a Lattice<T>,
    extent_iter: ExtentIterator,
}

impl<'a, T: Clone> Iterator for LatticeKeyValIterator<'a, T> {
    type Item = (Point, T);

    fn next(&mut self) -> Option<Self::Item> {
        self.extent_iter
            .next()
            .map(|p| (p, self.lattice.get_world(&p).clone()))
    }
}

pub struct ChunkedLatticeIterator<'a, T, I>
where
    I: Iterator<Item = &'a Lattice<T>>,
{
    full_extent: Extent,
    lattices: I,
    lattice_iter: Option<LatticeKeyValIterator<'a, T>>,
}

impl<'a, T, I> ChunkedLatticeIterator<'a, T, I>
where
    T: Clone,
    I: Iterator<Item = &'a Lattice<T>>,
{
    fn move_to_next_lattice(&mut self) {
        self.lattice_iter = self.lattices.next().map(|l| {
            let extent_iter = l.get_extent().intersection(&self.full_extent).into_iter();

            LatticeKeyValIterator {
                lattice: l,
                extent_iter,
            }
        });
    }
}

impl<'a, T, I> Iterator for ChunkedLatticeIterator<'a, T, I>
where
    T: Clone,
    I: Iterator<Item = &'a Lattice<T>>,
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
    fn test_periodic_indexer() {
        let extent = Extent::from_min_and_local_supremum([-1, -1, -1].into(), [3, 3, 3].into());
        let mut lattice = Lattice::fill_with_indexer(PeriodicYLevelsIndexer {}, extent, (0, 0, 0));
        for p in &extent {
            *lattice.get_mut_world(&p) = (p.x, p.y, p.z);
        }

        assert_eq!(*lattice.get_world(&[-1, -1, -1].into()), (-1, -1, -1));

        assert_eq!(*lattice.get_world(&[-2, -1, -1].into()), (1, -1, -1));
        assert_eq!(*lattice.get_world(&[-1, -2, -1].into()), (-1, 1, -1));
        assert_eq!(*lattice.get_world(&[-1, -1, -2].into()), (-1, -1, 1));

        assert_eq!(*lattice.get_world(&[-3, -1, -1].into()), (0, -1, -1));
        assert_eq!(*lattice.get_world(&[-1, -3, -1].into()), (-1, 0, -1));
        assert_eq!(*lattice.get_world(&[-1, -1, -3].into()), (-1, -1, 0));
    }

    #[test]
    fn test_local_point_from_index() {
        let sup = Point::new(10, 20, 30);
        let test_points = [
            Point::new(0, 0, 0),
            Point::new(1, 1, 1),
            Point::new(1, 2, 3),
            Point::new(4, 3, 2),
            Point::new(9, 9, 10),
            Point::new(9, 19, 29),
        ];

        for p in test_points.iter() {
            assert_eq!(
                *p,
                YLevelsIndexer::local_point_from_index(
                    &sup,
                    YLevelsIndexer::index_from_local_point(&sup, p)
                ),
            );
        }
    }

    #[test]
    fn test_serialized_extent_back_to_lattice() {
        let extent = Extent::from_min_and_local_supremum([-1, -1, -1].into(), [3, 3, 3].into());
        let mut lattice = Lattice::fill(extent, (0, 0, 0));
        for p in &extent {
            *lattice.get_mut_world(&p) = (p.x, p.y, p.z);
        }
        let orig_lattice = lattice.clone();

        let serial_extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 2, 2].into());
        let serialized = lattice.serialize_extent(&serial_extent);
        let serial_lattice = Lattice::new(serial_extent, serialized);
        Lattice::copy_extent(&serial_lattice, &mut lattice, &serial_extent);

        assert_eq!(orig_lattice, lattice);
    }

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

    #[test]
    fn test_apply_identity() {
        let matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let tfm = Transform { matrix };

        let extent = Extent::from_min_and_local_supremum([-1, -1, -1].into(), [3, 3, 3].into());
        let mut lattice = Lattice::fill(extent, (0, 0, 0));
        for p in &extent {
            *lattice.get_mut_world(&p) = (p.x, p.y, p.z);
        }

        let orig_lattice = lattice.clone();
        let tfm_lattice = lattice.apply_octahedral_transform(&tfm);

        assert_eq!(orig_lattice, tfm_lattice);
    }

    #[test]
    fn test_apply_x_reflection() {
        let matrix = [[-1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let tfm = Transform { matrix };

        let orig_extent =
            Extent::from_min_and_local_supremum([-2, -2, -2].into(), [3, 3, 3].into());
        let mut orig_lattice = Lattice::fill(orig_extent, (0, 0, 0));
        for p in &orig_extent {
            *orig_lattice.get_mut_world(&p) = (p.x, p.y, p.z);
        }

        let tfm_lattice = orig_lattice.apply_octahedral_transform(&tfm);

        let manual_tfm_extent =
            Extent::from_min_and_local_supremum([0, -2, -2].into(), [3, 3, 3].into());
        let mut manual_tfm_lattice = Lattice::fill(manual_tfm_extent, (0, 0, 0));
        for p in &manual_tfm_extent {
            *manual_tfm_lattice.get_mut_world(&p) = (-p.x, p.y, p.z);
        }

        assert_eq!(tfm_lattice, manual_tfm_lattice);
    }

    #[test]
    fn test_apply_rotation() {
        let matrix = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]];
        let tfm = Transform { matrix };

        let orig_extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), [1, 2, 3].into());
        let mut orig_lattice = Lattice::fill(orig_extent, (0, 0, 0));
        for p in &orig_extent {
            *orig_lattice.get_mut_world(&p) = (p.x, p.y, p.z);
        }

        let tfm_lattice = orig_lattice.apply_octahedral_transform(&tfm);

        let manual_tfm_extent =
            Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 1, 3].into());
        let mut manual_tfm_lattice = Lattice::fill(manual_tfm_extent, (0, 0, 0));
        for p in &manual_tfm_extent {
            // Have to do the inverse rotation here.
            *manual_tfm_lattice.get_mut_world(&p) = (-p.y, p.x, p.z);
        }

        assert_eq!(tfm_lattice, manual_tfm_lattice);
    }

    #[test]
    fn test_rotationally_symmetric_lattice_serializes_equivalently() {
        let orig_extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 2, 2].into());
        let bottom_extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 2, 1].into());
        let top_extent = Extent::from_min_and_local_supremum([0, 0, 1].into(), [2, 2, 1].into());
        let mut orig_lattice = Lattice::fill(orig_extent, 0);
        for p in &bottom_extent {
            *orig_lattice.get_mut_world(&p) = 1;
        }
        for p in &top_extent {
            *orig_lattice.get_mut_world(&p) = 2;
        }

        let orig_serial = orig_lattice.serialize_extent(&orig_extent);

        let matrix = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]];
        let tfm = Transform { matrix };

        // Apply successive rotations and make sure the serialization is always the same.
        let mut prev_lattice = orig_lattice;
        for _ in 0..4 {
            let tfm_lattice = prev_lattice.apply_octahedral_transform(&tfm);
            let tfm_serial = tfm_lattice.serialize_extent(&tfm_lattice.get_extent());
            assert_eq!(tfm_serial, orig_serial);
            prev_lattice = tfm_lattice;
        }
    }
}
