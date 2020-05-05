use crate::{Extent, Point};

use serde::{Deserialize, Serialize};

/// All of the possible symmetries of a Lattice (octahedron).
pub const OCTAHEDRAL_GROUP: [[[i32; 3]; 3]; 48] = [
    // Stationary X
    [[1, 0, 0], [0, 0, 1], [0, 1, 0]],
    [[1, 0, 0], [0, 0, -1], [0, 1, 0]],
    [[1, 0, 0], [0, 0, 1], [0, -1, 0]],
    [[1, 0, 0], [0, 0, -1], [0, -1, 0]],
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
    [[1, 0, 0], [0, 1, 0], [0, 0, -1]],
    [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
    // Mirror X
    [[-1, 0, 0], [0, 0, 1], [0, 1, 0]],
    [[-1, 0, 0], [0, 0, -1], [0, 1, 0]],
    [[-1, 0, 0], [0, 0, 1], [0, -1, 0]],
    [[-1, 0, 0], [0, 0, -1], [0, -1, 0]],
    [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
    [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
    // Stationary Y
    [[0, 0, 1], [0, 1, 0], [1, 0, 0]],
    [[0, 0, -1], [0, 1, 0], [1, 0, 0]],
    [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],
    [[0, 0, -1], [0, 1, 0], [-1, 0, 0]],
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[1, 0, 0], [0, 1, 0], [0, 0, -1]],
    [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
    // Mirror Y
    [[0, 0, 1], [0, -1, 0], [1, 0, 0]],
    [[0, 0, -1], [0, -1, 0], [1, 0, 0]],
    [[0, 0, 1], [0, -1, 0], [-1, 0, 0]],
    [[0, 0, -1], [0, -1, 0], [-1, 0, 0]],
    [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
    [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
    [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
    // Stationary Z
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
    [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
    [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
    [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
    [[0, -1, 0], [-1, 0, 0], [0, 0, 1]],
    // Mirror Z
    [[1, 0, 0], [0, 1, 0], [0, 0, -1]],
    [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
    [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
    [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
    [[0, 1, 0], [1, 0, 0], [0, 0, -1]],
    [[0, -1, 0], [1, 0, 0], [0, 0, -1]],
    [[0, 1, 0], [-1, 0, 0], [0, 0, -1]],
    [[0, -1, 0], [-1, 0, 0], [0, 0, -1]],
];

/// Useful for keeping voxels "upright."
pub const Z_STATIONARY_OCTAHEDRAL_GROUP: [[[i32; 3]; 3]; 8] = [
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
    [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
    [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
    [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
    [[0, -1, 0], [-1, 0, 0], [0, 0, 1]],
];

/// All octahedral without double mirrors.
pub const Z_STATIONARY_NO_DOUBLE_MIRROR_GROUP: [[[i32; 3]; 3]; 6] = [
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
    [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
    [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
    [[0, -1, 0], [-1, 0, 0], [0, 0, 1]],
];

/// Octahedral without rotations.
pub const Z_STATIONARY_MIRROR_GROUP: [[[i32; 3]; 3]; 4] = [
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
];

/// A linear map from (i32, i32, i32) to (i32, i32, i32).
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, Hash, PartialEq, Serialize)]
pub struct Transform {
    pub matrix: [[i32; 3]; 3],
}

impl Transform {
    pub fn apply_to_point(&self, p: &Point) -> Point {
        let Transform { matrix } = self;

        let x_map: Point = matrix[0].into();
        let y_map: Point = matrix[1].into();
        let z_map: Point = matrix[2].into();

        Point::new(x_map.dot(p), y_map.dot(p), z_map.dot(p))
    }

    pub fn apply_to_extent(&self, extent: &Extent) -> Extent {
        let corners = extent.get_world_corners();
        let tfm_corners: Vec<[i32; 3]> = corners
            .iter()
            .map(|c| self.apply_to_point(c).into())
            .collect();
        let tfm_min: Point = (*tfm_corners.iter().min().unwrap()).into();
        let tfm_max: Point = (*tfm_corners.iter().max().unwrap()).into();

        Extent::from_min_and_world_max(tfm_min, tfm_max)
    }

    pub fn is_octahedral(&self) -> bool {
        for matrix in OCTAHEDRAL_GROUP.iter() {
            if *matrix == self.matrix {
                return true;
            }
        }

        false
    }
}
