use crate::{Extent, ExtentIterator, Indexer, Point, StatelessIndexer, Transform, YLevelsIndexer};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Lattice<T, I = YLevelsIndexer> {
    indexer: I,
    extent: Extent,
    values: Vec<T>,
}

impl<T: Clone, I: Indexer> Lattice<T, I> {
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

    pub fn translate(&mut self, delta: &Point) {
        self.extent = self.extent + *delta;
    }

    pub fn fill_with_indexer(indexer: I, extent: Extent, init_val: T) -> Self {
        Lattice {
            extent,
            indexer,
            values: vec![init_val; extent.volume()],
        }
    }

    /// Returns a vec of the data in `extent`, ordered linearly by `I: Indexer`. A lattice
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

    pub fn copy_extent_to_position<S: From<T>, J: Indexer>(
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

    pub fn copy_extent<S: From<T>, J: Indexer>(
        src: &Self,
        dst: &mut Lattice<S, J>,
        extent: &Extent,
    ) {
        Self::copy_extent_to_position(src, dst, &extent.get_minimum(), extent)
    }

    pub fn copy_extent_into_new_lattice(&self, extent: &Extent) -> Self {
        let volume = extent.volume();
        let mut values = Vec::with_capacity(volume);
        unsafe { values.set_len(volume); }
        let mut copy = Lattice::new_with_indexer(*extent, self.indexer.clone(), values);
        Self::copy_extent(self, &mut copy, extent);

        copy
    }

    pub fn map_extent<S, F, J>(src: &Self, dst: &mut Lattice<S, J>, extent: &Extent, f: F)
    where
        F: Fn(&T) -> S,
        J: Indexer,
    {
        for p in extent {
            *dst.get_mut_world(&p) = f(src.get_world(&p));
        }
    }
}

impl<T, I: StatelessIndexer> Lattice<T, I> {
    pub fn new(extent: Extent, values: Vec<T>) -> Self {
        Lattice {
            extent,
            indexer: I::new(),
            values,
        }
    }

    pub fn new_at_origin(sup: Point, values: Vec<T>) -> Self {
        let extent = Extent::from_min_and_world_supremum([0, 0, 0].into(), sup);

        Self::new(extent, values)
    }
}

impl<T: Clone, I: StatelessIndexer> Lattice<T, I> {
    pub fn fill(extent: Extent, init_val: T) -> Self {
        Self::fill_with_indexer(I::new(), extent, init_val)
    }
}

impl<T, I: Indexer> Lattice<T, I> {
    pub fn new_with_indexer(extent: Extent, indexer: I, values: Vec<T>) -> Self {
        Lattice {
            extent,
            indexer,
            values,
        }
    }

    pub fn get_indexer(&self) -> &I {
        &self.indexer
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

pub struct LatticeKeyValIterator<'a, T> {
    lattice: &'a Lattice<T>,
    extent_iter: ExtentIterator,
}

impl<'a, T> LatticeKeyValIterator<'a, T> {
    pub fn new(lattice: &'a Lattice<T>, extent_iter: ExtentIterator) -> Self {
        Self { lattice, extent_iter }
    }
}

impl<'a, T: Clone> Iterator for LatticeKeyValIterator<'a, T> {
    type Item = (Point, T);

    fn next(&mut self) -> Option<Self::Item> {
        self.extent_iter
            .next()
            .map(|p| (p, self.lattice.get_world(&p).clone()))
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
    use crate::PeriodicYLevelsIndexer;

    #[test]
    fn test_periodic_indexer() {
        let extent = Extent::from_min_and_local_supremum([-1, -1, -1].into(), [3, 3, 3].into());
        let mut lattice = Lattice::<_, PeriodicYLevelsIndexer>::fill(extent, (0, 0, 0));
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
        let mut lattice = Lattice::<_, YLevelsIndexer>::fill(extent, (0, 0, 0));
        for p in &extent {
            *lattice.get_mut_world(&p) = (p.x, p.y, p.z);
        }
        let orig_lattice = lattice.clone();

        let serial_extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 2, 2].into());
        let serialized = lattice.serialize_extent(&serial_extent);
        let serial_lattice = Lattice::<_, YLevelsIndexer>::new(serial_extent, serialized);
        Lattice::copy_extent(&serial_lattice, &mut lattice, &serial_extent);

        assert_eq!(orig_lattice, lattice);
    }

    #[test]
    fn test_apply_identity() {
        let matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let tfm = Transform { matrix };

        let extent = Extent::from_min_and_local_supremum([-1, -1, -1].into(), [3, 3, 3].into());
        let mut lattice = Lattice::<_, YLevelsIndexer>::fill(extent, (0, 0, 0));
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
        let mut orig_lattice = Lattice::<_, YLevelsIndexer>::fill(orig_extent, (0, 0, 0));
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
        let mut orig_lattice = Lattice::<_, YLevelsIndexer>::fill(orig_extent, (0, 0, 0));
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
        let mut orig_lattice = Lattice::<_, YLevelsIndexer>::fill(orig_extent, 0);
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
