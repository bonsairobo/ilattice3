use crate::{
    copy_extent, lattice_map::LatticeMapKeyValIterator, prelude::*, Extent, Indexer,
    StatelessIndexer, Transform, YLevelsIndexer,
};

use serde::{Deserialize, Serialize};

/// A map from points in an extent to some kind of data `T`.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct VecLatticeMap<T, I = YLevelsIndexer> {
    // Works if I: Default, which is true for I: StatelessIndexer.
    #[serde(skip_deserializing)]
    #[serde(skip_serializing)]
    indexer: I,

    extent: Extent,
    values: Vec<T>,
}

impl<T, I> GetExtent for VecLatticeMap<T, I> {
    fn get_extent(&self) -> &Extent {
        &self.extent
    }
}

impl<T, I: Indexer> GetLocal<T> for VecLatticeMap<T, I>
where
    T: Clone,
{
    fn get_local(&self, p: &Point) -> T {
        self.get_linear_ref(self.index_from_local_point(p)).clone()
    }
}

impl<T, I: Indexer> GetLocalRef<T> for VecLatticeMap<T, I> {
    fn get_local_ref(&self, p: &Point) -> &T {
        self.get_linear_ref(self.index_from_local_point(p))
    }
}

impl<T, I: Indexer> GetLocalRefMut<T> for VecLatticeMap<T, I> {
    fn get_local_ref_mut(&mut self, p: &Point) -> &mut T {
        self.get_linear_ref_mut(self.index_from_local_point(p))
    }
}

impl<T, I: Indexer> VecLatticeMap<T, I> {
    pub fn new_with_indexer(extent: Extent, indexer: I, values: Vec<T>) -> Self {
        VecLatticeMap {
            extent,
            indexer,
            values,
        }
    }

    pub fn get_indexer(&self) -> &I {
        &self.indexer
    }

    pub fn get_linear_ref(&self, index: usize) -> &T {
        &self.values[index]
    }

    pub fn get_linear_ref_mut(&mut self, index: usize) -> &mut T {
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

    pub fn index_from_world_point(&self, p: &Point) -> usize {
        self.index_from_local_point(&self.extent.local_point_from_world_point(p))
    }

    pub fn translate(&mut self, delta: &Point) {
        self.extent = self.extent + *delta;
    }

    pub fn set_minimum(&mut self, new_min: &Point) {
        self.extent = self.extent.with_minimum(*new_min);
    }

    pub fn map<F, S>(&self, f: F) -> VecLatticeMap<S, I>
    where
        F: Fn(&T) -> S,
    {
        VecLatticeMap::new_with_indexer(
            *self.get_extent(),
            self.indexer.clone(),
            self.values.iter().map(f).collect(),
        )
    }

    pub fn iter(&self) -> LatticeMapKeyValIterator<VecLatticeMap<T, I>, T> {
        LatticeMapKeyValIterator::new(self, self.extent.into_iter())
    }
}

impl<T, I: StatelessIndexer> VecLatticeMap<T, I> {
    pub fn new(extent: Extent, values: Vec<T>) -> Self {
        VecLatticeMap {
            extent,
            indexer: I::default(),
            values,
        }
    }

    pub fn new_at_origin(sup: Point, values: Vec<T>) -> Self {
        let extent = Extent::from_min_and_world_supremum([0, 0, 0].into(), sup);

        Self::new(extent, values)
    }
}

impl<T: Clone, I: Indexer> VecLatticeMap<T, I> {
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
        let mut tfm_map = Self::new_with_indexer(tfm_extent, self.indexer.clone(), new_values);

        // PERF: this is not the most efficient, but it is very simple.
        for p in extent {
            let tfm_p = tfm.apply_to_point(&p);
            *tfm_map.get_world_ref_mut(&tfm_p) = self.get_world(&p);
        }

        tfm_map
    }

    pub fn fill_with_indexer(indexer: I, extent: Extent, init_val: T) -> Self {
        VecLatticeMap {
            extent,
            indexer,
            values: vec![init_val; extent.volume()],
        }
    }

    /// Returns a vec of the data in `extent`, ordered linearly by `I: Indexer`. A map
    /// can be recreated from the vec using `VecLatticeMap::<T, I>::new_with_indexer`.
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
            data[i] = self.get_world(&p);
        }

        data
    }

    pub fn copy_extent_into_new_map(&self, extent: &Extent) -> Self {
        let volume = extent.volume();
        let mut values = Vec::with_capacity(volume);
        unsafe {
            values.set_len(volume);
        }
        let mut copy = VecLatticeMap::new_with_indexer(*extent, self.indexer.clone(), values);
        copy_extent(self, &mut copy, extent);

        copy
    }
}

impl<T: Clone, I: StatelessIndexer> VecLatticeMap<T, I> {
    pub fn fill(extent: Extent, init_val: T) -> Self {
        Self::fill_with_indexer(I::default(), extent, init_val)
    }

    pub fn copy_from_map<V>(map: &V, extent: &Extent) -> Self
    where
        V: GetWorld<T>,
    {
        let volume = extent.volume();
        let mut values = Vec::with_capacity(volume);
        unsafe {
            values.set_len(volume);
        }
        let mut copy = VecLatticeMap::new(*extent, values);
        copy_extent(map, &mut copy, extent);

        copy
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
        let mut map = VecLatticeMap::<_, PeriodicYLevelsIndexer>::fill(extent, (0, 0, 0));
        for p in &extent {
            *map.get_world_ref_mut(&p) = (p.x, p.y, p.z);
        }

        assert_eq!(map.get_world(&[-1, -1, -1].into()), (-1, -1, -1));

        assert_eq!(map.get_world(&[-2, -1, -1].into()), (1, -1, -1));
        assert_eq!(map.get_world(&[-1, -2, -1].into()), (-1, 1, -1));
        assert_eq!(map.get_world(&[-1, -1, -2].into()), (-1, -1, 1));

        assert_eq!(map.get_world(&[-3, -1, -1].into()), (0, -1, -1));
        assert_eq!(map.get_world(&[-1, -3, -1].into()), (-1, 0, -1));
        assert_eq!(map.get_world(&[-1, -1, -3].into()), (-1, -1, 0));
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
    fn test_serialized_extent_back_to_map() {
        let extent = Extent::from_min_and_local_supremum([-1, -1, -1].into(), [3, 3, 3].into());
        let mut map = VecLatticeMap::<_, YLevelsIndexer>::fill(extent, (0, 0, 0));
        for p in &extent {
            *map.get_world_ref_mut(&p) = (p.x, p.y, p.z);
        }
        let orig_map = map.clone();

        let serial_extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 2, 2].into());
        let serialized = map.serialize_extent(&serial_extent);
        let serial_map = VecLatticeMap::<_, YLevelsIndexer>::new(serial_extent, serialized);
        copy_extent(&serial_map, &mut map, &serial_extent);

        assert_eq!(orig_map, map);
    }

    #[test]
    fn test_apply_identity() {
        let matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let tfm = Transform { matrix };

        let extent = Extent::from_min_and_local_supremum([-1, -1, -1].into(), [3, 3, 3].into());
        let mut map = VecLatticeMap::<_, YLevelsIndexer>::fill(extent, (0, 0, 0));
        for p in &extent {
            *map.get_world_ref_mut(&p) = (p.x, p.y, p.z);
        }

        let orig_map = map.clone();
        let tfm_map = map.apply_octahedral_transform(&tfm);

        assert_eq!(orig_map, tfm_map);
    }

    #[test]
    fn test_apply_x_reflection() {
        let matrix = [[-1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let tfm = Transform { matrix };

        let orig_extent =
            Extent::from_min_and_local_supremum([-2, -2, -2].into(), [3, 3, 3].into());
        let mut orig_map = VecLatticeMap::<_, YLevelsIndexer>::fill(orig_extent, (0, 0, 0));
        for p in &orig_extent {
            *orig_map.get_world_ref_mut(&p) = (p.x, p.y, p.z);
        }

        let tfm_map = orig_map.apply_octahedral_transform(&tfm);

        let manual_tfm_extent =
            Extent::from_min_and_local_supremum([0, -2, -2].into(), [3, 3, 3].into());
        let mut manual_tfm_map = VecLatticeMap::fill(manual_tfm_extent, (0, 0, 0));
        for p in &manual_tfm_extent {
            *manual_tfm_map.get_world_ref_mut(&p) = (-p.x, p.y, p.z);
        }

        assert_eq!(tfm_map, manual_tfm_map);
    }

    #[test]
    fn test_apply_rotation() {
        let matrix = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]];
        let tfm = Transform { matrix };

        let orig_extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), [1, 2, 3].into());
        let mut orig_map = VecLatticeMap::<_, YLevelsIndexer>::fill(orig_extent, (0, 0, 0));
        for p in &orig_extent {
            *orig_map.get_world_ref_mut(&p) = (p.x, p.y, p.z);
        }

        let tfm_map = orig_map.apply_octahedral_transform(&tfm);

        let manual_tfm_extent =
            Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 1, 3].into());
        let mut manual_tfm_map = VecLatticeMap::fill(manual_tfm_extent, (0, 0, 0));
        for p in &manual_tfm_extent {
            // Have to do the inverse rotation here.
            *manual_tfm_map.get_world_ref_mut(&p) = (-p.y, p.x, p.z);
        }

        assert_eq!(tfm_map, manual_tfm_map);
    }

    #[test]
    fn test_rotationally_symmetric_map_serializes_equivalently() {
        let orig_extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 2, 2].into());
        let bottom_extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 2, 1].into());
        let top_extent = Extent::from_min_and_local_supremum([0, 0, 1].into(), [2, 2, 1].into());
        let mut orig_map = VecLatticeMap::<_, YLevelsIndexer>::fill(orig_extent, 0);
        for p in &bottom_extent {
            *orig_map.get_world_ref_mut(&p) = 1;
        }
        for p in &top_extent {
            *orig_map.get_world_ref_mut(&p) = 2;
        }

        let orig_serial = orig_map.serialize_extent(&orig_extent);

        let matrix = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]];
        let tfm = Transform { matrix };

        // Apply successive rotations and make sure the serialization is always the same.
        let mut prev_map = orig_map;
        for _ in 0..4 {
            let tfm_map = prev_map.apply_octahedral_transform(&tfm);
            let tfm_serial = tfm_map.serialize_extent(&tfm_map.get_extent());
            assert_eq!(tfm_serial, orig_serial);
            prev_map = tfm_map;
        }
    }
}
