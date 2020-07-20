use crate::{Extent, GetLinearRef, HasIndexer, Indexer, IsEmpty, Point};

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
    for (i, offset) in FACE_ADJACENT_OFFSETS[0..3].iter().enumerate() {
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
            let is_empty = voxels.get_linear_ref(p_linear_offset as usize).is_empty();
            if is_empty {
                surface_points.push(p + min);
                break;
            }
        }
    }

    surface_points
}

const FACE_ADJACENT_OFFSETS: [Point; 6] = [
    Point { x: 1, y: 0, z: 0 },
    Point { x: 0, y: 1, z: 0 },
    Point { x: 0, y: 0, z: 1 },
    Point { x: -1, y: 0, z: 0 },
    Point { x: 0, y: -1, z: 0 },
    Point { x: 0, y: 0, z: -1 },
];
