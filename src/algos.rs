use crate::{Extent, GetLinearRef, HasIndexer, Indexer, IsEmpty, Point, FACE_ADJACENT_OFFSETS};

/// Returns the voxels that are considered "visible," i.e. they are solid and face-adjacent to an
/// empty voxel. Voxels on the boundary of `extent` will not be considered candidates for being
/// "visible," but they will be used when checking adjacent voxels.
pub fn find_surface_voxels<V, T, I>(voxels: &V, extent: &Extent) -> Vec<Point>
where
    V: GetLinearRef<Data = T> + HasIndexer<Indexer = I>,
    T: IsEmpty,
    I: Indexer,
{
    let min = extent.get_minimum();
    let sup = extent.get_local_supremum();

    // Precompute the offsets for adjacency checks. Some of these usize values will be considered
    // "negative," although they have wrapped around.
    let mut linear_offsets = [0; 6];
    for (i, offset) in FACE_ADJACENT_OFFSETS.iter().enumerate() {
        linear_offsets[i] = I::index_from_local_point(sup, offset);
    }

    let mut surface_points = Vec::new();
    let iter_extent = extent.with_minimum([0, 0, 0].into()).radial_grow(-1);
    for p in iter_extent {
        let p_linear = I::index_from_local_point(sup, &p);
        if voxels.get_linear_ref(p_linear).is_empty() {
            continue;
        }
        for linear_offset in linear_offsets.iter() {
            let p_linear_offset = p_linear.wrapping_add(*linear_offset);
            let is_empty = voxels.get_linear_ref(p_linear_offset).is_empty();
            if is_empty {
                surface_points.push(p + min);
                break;
            }
        }
    }

    surface_points
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        fill_extent, prelude::*, test_util::assert_elements_eq, VecLatticeMap, YLevelsIndexer,
    };

    #[derive(Clone)]
    struct Voxel(bool);

    impl IsEmpty for Voxel {
        fn is_empty(&self) -> bool {
            !self.0
        }
    }

    #[test]
    fn find_surface_voxels_cube_side_length_3() {
        let mut voxels: VecLatticeMap<_, YLevelsIndexer> = VecLatticeMap::fill(
            Extent::from_min_and_local_supremum([0, 0, 0].into(), [5, 5, 5].into()),
            Voxel(false),
        );

        let center = [2, 2, 2].into();
        let solid_extent = Extent::from_center_and_radius(center, 1);
        fill_extent(&mut voxels, &solid_extent, Voxel(true));

        let surface_points = find_surface_voxels(&voxels, voxels.get_extent());

        // Should exclude the center point.
        let expected_surface_points = solid_extent.into_iter().filter(|p| *p != center).collect();
        assert_elements_eq(&surface_points, &expected_surface_points);
    }
}
