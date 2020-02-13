use num::Integer;
use serde::{Deserialize, Serialize};
use std::cmp::{max, min, Ordering};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

// TODO: SIMD?

#[derive(Copy, Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
pub struct Point {
    pub x: i32,
    pub y: i32,
    pub z: i32,
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}, {}, {})", self.x, self.y, self.z)
    }
}

impl From<Point> for [f32; 3] {
    fn from(p: Point) -> [f32; 3] {
        [p.x as f32, p.y as f32, p.z as f32]
    }
}

impl From<Point> for [i32; 3] {
    fn from(p: Point) -> [i32; 3] {
        [p.x, p.y, p.z]
    }
}

impl From<[i32; 3]> for Point {
    fn from(other: [i32; 3]) -> Self {
        Point::new(other[0], other[1], other[2])
    }
}

impl From<(i32, i32, i32)> for Point {
    fn from(other: (i32, i32, i32)) -> Self {
        let (x, y, z) = other;

        Point::new(x, y, z)
    }
}

impl Point {
    pub fn new(x: i32, y: i32, z: i32) -> Self {
        Self { x, y, z }
    }

    pub fn zero() -> Self {
        [0, 0, 0].into()
    }

    pub fn xyz(&self) -> [i32; 3] {
        [self.x, self.y, self.z]
    }

    pub fn zxy(&self) -> [i32; 3] {
        [self.z, self.x, self.y]
    }

    pub fn yzx(&self) -> [i32; 3] {
        [self.y, self.z, self.x]
    }

    pub fn xzy(&self) -> [i32; 3] {
        [self.x, self.z, self.y]
    }

    pub fn yxz(&self) -> [i32; 3] {
        [self.y, self.x, self.z]
    }

    pub fn zyx(&self) -> [i32; 3] {
        [self.z, self.y, self.x]
    }

    pub fn xy(&self) -> [i32; 2] {
        [self.x, self.y]
    }

    pub fn zx(&self) -> [i32; 2] {
        [self.z, self.x]
    }

    pub fn yz(&self) -> [i32; 2] {
        [self.y, self.z]
    }

    pub fn xz(&self) -> [i32; 2] {
        [self.x, self.z]
    }

    pub fn yx(&self) -> [i32; 2] {
        [self.y, self.x]
    }

    pub fn zy(&self) -> [i32; 2] {
        [self.z, self.y]
    }

    pub fn map_components(&self, f: &impl Fn(i32) -> i32) -> Self {
        [f(self.x), f(self.y), f(self.z)].into()
    }

    pub fn dot(&self, other: &Self) -> i32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn join(&self, other: &Self) -> Self {
        [
            max(self.x, other.x),
            max(self.y, other.y),
            max(self.z, other.z),
        ]
        .into()
    }

    pub fn meet(&self, other: &Self) -> Self {
        [
            min(self.x, other.x),
            min(self.y, other.y),
            min(self.z, other.z),
        ]
        .into()
    }

    pub fn norm(&self) -> f32 {
        (self.dot(self) as f32).sqrt()
    }
}

/// This particular partial order allows us to say that a lattice L contains a point iff
/// p is GEQ the minimum of L and p is LEQ the maximum of L.
impl PartialOrd for Point {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self < other {
            Some(Ordering::Less)
        } else if self > other {
            Some(Ordering::Greater)
        } else if self.x == other.x && self.y == other.y && self.z == other.z {
            Some(Ordering::Equal)
        } else {
            None
        }
    }

    fn lt(&self, other: &Self) -> bool {
        self.x < other.x && self.y < other.y && self.z < other.z
    }

    fn gt(&self, other: &Self) -> bool {
        self.x > other.x && self.y > other.y && self.z > other.z
    }

    fn le(&self, other: &Self) -> bool {
        self.x <= other.x && self.y <= other.y && self.z <= other.z
    }

    fn ge(&self, other: &Self) -> bool {
        self.x >= other.x && self.y >= other.y && self.z >= other.z
    }
}

impl Add for Point {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        [self.x + other.x, self.y + other.y, self.z + other.z].into()
    }
}

impl Neg for Point {
    type Output = Self;

    fn neg(self) -> Self {
        [-self.x, -self.y, -self.z].into()
    }
}

impl Sub for Point {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}

impl Mul<i32> for Point {
    type Output = Self;

    fn mul(self, rhs: i32) -> Self {
        [rhs * self.x, rhs * self.y, rhs * self.z].into()
    }
}

impl Mul<Point> for Point {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        [other.x * self.x, other.y * self.y, other.z * self.z].into()
    }
}

impl Div<i32> for Point {
    type Output = Self;

    fn div(self, rhs: i32) -> Self {
        [
            self.x.div_floor(&rhs),
            self.y.div_floor(&rhs),
            self.z.div_floor(&rhs),
        ]
        .into()
    }
}

impl Div<Point> for Point {
    type Output = Self;

    fn div(self, rhs: Point) -> Self {
        [
            self.x.div_floor(&rhs.x),
            self.y.div_floor(&rhs.y),
            self.z.div_floor(&rhs.z),
        ]
        .into()
    }
}
