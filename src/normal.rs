use crate::Point;

use enum_primitive_derive::Primitive;
use num_traits::cast::FromPrimitive;

/// The integer order of directions is often used for indexing and arithmetic direction
/// transformations.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, Primitive)]
pub enum Direction {
    NegX = 0,
    PosX = 1,
    NegY = 2,
    PosY = 3,
    NegZ = 4,
    PosZ = 5,
}

impl Direction {
    pub fn is_negative(self) -> bool {
        self as u32 % 2 == 0
    }

    pub fn is_x(self) -> bool {
        self as u32 / 2 == 0
    }

    pub fn is_y(self) -> bool {
        self as u32 / 2 == 1
    }

    pub fn is_z(self) -> bool {
        self as u32 / 2 == 2
    }

    pub fn negate(self) -> Self {
        let n = self as u32;
        let div2 = n / 2;
        let mod2 = n % 2;

        Direction::from_u32(2 * div2 + (mod2 + 1) % 2).expect("Bad direction index")
    }

    pub fn positive(self) -> Self {
        if self.is_negative() {
            self.negate()
        } else {
            self
        }
    }
}

/// Allows for indexing by `Direction`.
#[derive(Clone, Copy)]
pub struct DirectionIndex<T> {
    /// Be careful to make each the value corresponds to the correct `Direction`.
    pub values: [T; 6],
}

impl<T> DirectionIndex<T> {
    pub fn new(values: [T; 6]) -> Self {
        DirectionIndex { values }
    }

    pub fn get(&self, direction: Direction) -> &T {
        &self.values[direction as usize]
    }

    pub fn get_mut(&mut self, direction: Direction) -> &mut T {
        &mut self.values[direction as usize]
    }

    pub fn iter(&self) -> impl Iterator<Item = (Direction, &T)> {
        self.values.iter().enumerate().map(|(i, v)| {
            (
                Direction::from_usize(i).expect("Bad index for direction"),
                v,
            )
        })
    }
}

impl DirectionIndex<i32> {
    pub fn zeroes() -> Self {
        DirectionIndex { values: [0; 6] }
    }

    pub fn vector_for_direction(&self, direction: Direction) -> Point {
        ALL_NORMALS[direction as usize] * *self.get(direction)
    }

    pub fn min_vector(&self) -> (Point, Direction) {
        let (i, _min_dist) = self
            .values
            .iter()
            .enumerate()
            .min_by_key(|&(_, d)| d)
            .expect("Impossibru!");
        let direction = Direction::from_usize(i).expect("Bad direction index");

        (self.vector_for_direction(direction), direction)
    }

    pub fn positive_point(&self) -> Point {
        [
            *self.get(Direction::PosX),
            *self.get(Direction::PosY),
            *self.get(Direction::PosZ),
        ]
        .into()
    }

    pub fn negative_point(&self) -> Point {
        [
            -*self.get(Direction::NegX),
            -*self.get(Direction::NegY),
            -*self.get(Direction::NegZ),
        ]
        .into()
    }
}

impl From<Point> for Direction {
    fn from(n: Point) -> Self {
        debug_assert!(
            (n.x != 0 && n.y == 0 && n.z == 0)
                || (n.x == 0 && n.y != 0 && n.z == 0)
                || (n.x == 0 && n.y == 0 && n.z != 0)
        );

        if n.x < 0 {
            Direction::NegX
        } else if n.x > 0 {
            Direction::PosX
        } else if n.y < 0 {
            Direction::NegY
        } else if n.y > 0 {
            Direction::PosY
        } else if n.z < 0 {
            Direction::NegZ
        } else if n.z > 0 {
            Direction::PosZ
        } else {
            panic!("Bad normal vector")
        }
    }
}

pub const ALL_DIRECTIONS: [Direction; 6] = [
    Direction::NegX,
    Direction::PosX,
    Direction::NegY,
    Direction::PosY,
    Direction::NegZ,
    Direction::PosZ,
];

pub const ALL_NORMALS: [Point; 6] = [
    Point { x: -1, y: 0, z: 0 },
    Point { x: 1, y: 0, z: 0 },
    Point { x: 0, y: -1, z: 0 },
    Point { x: 0, y: 1, z: 0 },
    Point { x: 0, y: 0, z: -1 },
    Point { x: 0, y: 0, z: 1 },
];

/// Returns vector pointing in direction.
impl From<Direction> for Point {
    fn from(d: Direction) -> Self {
        ALL_NORMALS[d as usize]
    }
}

/// Helps convert between the two representations of axis directions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Normal {
    Vector(Point),
    Axis(Direction),
}

impl From<Normal> for Point {
    fn from(other: Normal) -> Self {
        match other {
            Normal::Vector(p) => p,
            Normal::Axis(a) => a.into(),
        }
    }
}

impl From<Normal> for Direction {
    fn from(other: Normal) -> Self {
        match other {
            Normal::Vector(v) => v.into(),
            Normal::Axis(a) => a,
        }
    }
}

/// Unit vectors spanning the plane of the quad.
pub struct PlaneSpanInfo {
    pub u: Point,
    pub v: Point,
}

impl Normal {
    pub fn as_axis(self) -> Self {
        Normal::Axis(self.into())
    }

    pub fn as_vector(self) -> Self {
        Normal::Vector(self.into())
    }

    /// Returns positive unit vectors spanning the plane designated by this normal.
    pub fn get_plane_span_info(&self) -> PlaneSpanInfo {
        let mut n = Point::from(*self);
        n = n.map_components(&|a| a.abs());

        PlaneSpanInfo {
            u: n.zxy().into(),
            v: n.yzx().into(),
        }
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

    #[test]
    fn test_direction_from_point() {
        let all_dirs: Vec<_> = ALL_DIRECTIONS
            .iter()
            .map(|d| Point::from(*d))
            .map(|p| Direction::from(p))
            .collect();

        assert_eq!(all_dirs, ALL_DIRECTIONS);
    }

    #[test]
    fn test_negative_directions() {
        assert_eq!(Direction::NegX.negate(), Direction::PosX);
        assert_eq!(Direction::PosX.negate(), Direction::NegX);
        assert_eq!(Direction::NegY.negate(), Direction::PosY);
        assert_eq!(Direction::PosY.negate(), Direction::NegY);
        assert_eq!(Direction::NegZ.negate(), Direction::PosZ);
        assert_eq!(Direction::PosZ.negate(), Direction::NegZ);
    }
}
