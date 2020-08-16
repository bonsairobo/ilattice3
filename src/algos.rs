use crate::{point::FACE_ADJACENT_OFFSETS, prelude::*};

/// Find all points inside of `extent` that are path-connected to `seed` and satisfy `predicate`.
/// `extent` and `seed` should be given in the world coordinates of `voxels`. Returns points in
/// world coordinates.
pub fn flood_fill_world<V, T, I>(
    voxels: &V,
    seed: Point,
    extent: &Extent,
    predicate: impl Fn(&T) -> bool,
) -> Vec<Point>
where
    V: GetExtent + GetLinearRef<Data = T> + HasIndexer<Indexer = I>,
    I: Indexer,
{
    let voxels_min = voxels.get_extent().get_minimum();
    let local_extent = *extent - voxels_min;
    let local_seed = seed - voxels_min;

    let local_filled = flood_fill_local(voxels, local_seed, &local_extent, predicate);
    let world_filled = local_filled
        .into_iter()
        .map(|p| p + extent.get_minimum())
        .collect();

    world_filled
}

// PERF: Could be faster with scanlines.
/// Find all points inside of `local_extent` that are path-connected to `local_seed` and satisfy
/// `predicate`. `local_extent` and `local_seed` should be given in the local coordinates of
/// `voxels`. Returns points in local coordinates.
pub fn flood_fill_local<V, T, I>(
    voxels: &V,
    local_seed: Point,
    local_extent: &Extent,
    predicate: impl Fn(&T) -> bool,
) -> Vec<Point>
where
    V: GetExtent + GetLinearRef<Data = T> + HasIndexer<Indexer = I>,
    I: Indexer,
{
    let mut found = Vec::new();

    let voxels_sup = voxels.get_extent().get_local_supremum();

    // Precompute the offsets for adjacency checks. Some of these usize values will be considered
    // "negative," although they have wrapped around.
    let mut linear_offsets = [0; 6];
    I::linear_strides(voxels_sup, &FACE_ADJACENT_OFFSETS, &mut linear_offsets);

    let min = local_extent.get_minimum();
    let linear_min = I::index_from_local_point(&voxels_sup, &min);

    let linear_seed = I::index_from_local_point(voxels_sup, &local_seed);

    let mut visited = vec![false; local_extent.volume()];
    let mut stack = vec![(local_seed, linear_seed)];

    while let Some((p, p_linear)) = stack.pop() {
        let visited_index = p_linear - linear_min;
        if visited[visited_index] {
            continue;
        }
        visited[visited_index] = true;
        if !predicate(voxels.get_linear_ref(p_linear)) {
            continue;
        }
        found.push(p);

        for (offset, linear_offset) in FACE_ADJACENT_OFFSETS.iter().zip(linear_offsets.iter()) {
            let p_offset = p + *offset;
            if local_extent.contains_world(&p_offset) {
                stack.push((p_offset, p_linear.wrapping_add(*linear_offset)))
            }
        }
    }

    found
}

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
    I::linear_strides(sup, &FACE_ADJACENT_OFFSETS, &mut linear_offsets);

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

// ████████╗███████╗███████╗████████╗███████╗
// ╚══██╔══╝██╔════╝██╔════╝╚══██╔══╝██╔════╝
//    ██║   █████╗  ███████╗   ██║   ███████╗
//    ██║   ██╔══╝  ╚════██║   ██║   ╚════██║
//    ██║   ███████╗███████║   ██║   ███████║
//    ╚═╝   ╚══════╝╚══════╝   ╚═╝   ╚══════╝

#[cfg(test)]
mod test {
    use super::*;
    use crate::{fill_extent, test_util::assert_elements_eq, VecLatticeMap, YLevelsIndexer};

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

    #[test]
    fn flood_fill_around_cube() {
        let extent = Extent::from_min_and_local_supremum([1, 1, 1].into(), [5, 5, 5].into());
        let mut voxels: VecLatticeMap<_, YLevelsIndexer> =
            VecLatticeMap::fill(extent, Voxel(false));

        let cube_center = [3, 3, 3].into();
        let cube_extent = Extent::from_center_and_radius(cube_center, 1);
        fill_extent(&mut voxels, &cube_extent, Voxel(true));

        let seed = [1, 1, 1].into();
        let predicate = |val: &Voxel| !val.0;
        let filled_points = flood_fill_world(&voxels, seed, voxels.get_extent(), &predicate);

        // Should exclude the center point.
        let expected_filled_points = extent.iter_boundary_points().collect();
        assert_elements_eq(&filled_points, &expected_filled_points);
    }
}
