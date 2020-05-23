use crate::{Extent, Point};

use serde::{Deserialize, Serialize};

pub type Matrix = [[i32; 3]; 3];

/// All of the possible symmetries of a Lattice (octahedron).
pub const OCTAHEDRAL_GROUP: [Matrix; 48] = [
    // (x, y, z)
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
    [[1, 0, 0], [0, 1, 0], [0, 0, -1]],
    [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
    [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
    [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],

    // (x, z, y)
    [[1, 0, 0], [0, 0, 1], [0, 1, 0]],
    [[1, 0, 0], [0, 0, -1], [0, 1, 0]],
    [[1, 0, 0], [0, 0, 1], [0, -1, 0]],
    [[1, 0, 0], [0, 0, -1], [0, -1, 0]],
    [[-1, 0, 0], [0, 0, 1], [0, 1, 0]],
    [[-1, 0, 0], [0, 0, 1], [0, 1, 0]],
    [[-1, 0, 0], [0, 0, 1], [0, 1, 0]],
    [[-1, 0, 0], [0, 0, 1], [0, 1, 0]],

    // (y, x, z)
    [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
    [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
    [[0, 1, 0], [1, 0, 0], [0, 0, -1]],
    [[0, 1, 0], [-1, 0, 0], [0, 0, -1]],
    [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
    [[0, -1, 0], [-1, 0, 0], [0, 0, 1]],
    [[0, -1, 0], [1, 0, 0], [0, 0, -1]],
    [[0, -1, 0], [-1, 0, 0], [0, 0, -1]],

    // (y, z, x)
    [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
    [[0, 1, 0], [0, 0, -1], [1, 0, 0]],
    [[0, 1, 0], [0, 0, 1], [-1, 0, 0]],
    [[0, 1, 0], [0, 0, -1], [-1, 0, 0]],
    [[0, -1, 0], [0, 0, 1], [1, 0, 0]],
    [[0, -1, 0], [0, 0, -1], [1, 0, 0]],
    [[0, -1, 0], [0, 0, 1], [-1, 0, 0]],
    [[0, -1, 0], [0, 0, -1], [-1, 0, 0]],

    // (z, y, x)
    [[0, 0, 1], [0, 1, 0], [1, 0, 0]],
    [[0, 0, 1], [0, -1, 0], [1, 0, 0]],
    [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],
    [[0, 0, 1], [0, -1, 0], [-1, 0, 0]],
    [[0, 0, -1], [0, 1, 0], [1, 0, 0]],
    [[0, 0, -1], [0, -1, 0], [1, 0, 0]],
    [[0, 0, -1], [0, 1, 0], [-1, 0, 0]],
    [[0, 0, -1], [0, -1, 0], [-1, 0, 0]],

    // (z, x, y)
    [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
    [[0, 0, 1], [-1, 0, 0], [0, 1, 0]],
    [[0, 0, 1], [1, 0, 0], [0, -1, 0]],
    [[0, 0, 1], [-1, 0, 0], [0, -1, 0]],
    [[0, 0, -1], [1, 0, 0], [0, 1, 0]],
    [[0, 0, -1], [-1, 0, 0], [0, 1, 0]],
    [[0, 0, -1], [1, 0, 0], [0, -1, 0]],
    [[0, 0, -1], [-1, 0, 0], [0, -1, 0]],
];

/// Useful for keeping voxels "upright."
pub const Z_STATIONARY_OCTAHEDRAL_GROUP: [Matrix; 8] = [
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
    [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
    [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
    [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
    [[0, -1, 0], [-1, 0, 0], [0, 0, 1]],
];

/// A linear map from (i32, i32, i32) to (i32, i32, i32).
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, Hash, PartialEq, Serialize)]
pub struct Transform {
    pub matrix: Matrix,
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
