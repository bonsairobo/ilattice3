use crate::{DirectionIndex, Point};

use std::ops::{Add, Sub};

/// A Cartesian product of 3 integer ranges: `[x_min..x_max] * [y_min..y_max] * [z_min..z_max]`.
#[derive(Copy, Clone, Debug, Eq, Hash, PartialEq)]
pub struct Extent {
    /// Point in the extent that's lesser than all other points in the extent. In world coordinates,
    /// since the local minimum is always (0, 0, 0).
    minimum: Point,

    /// A strict supremum is the least point that's strictly greater than all points in the extent.
    world_sup: Point,

    /// By "local," we mean relative to the minimum, so the local strict supremum effectively
    /// encodes the size of the extent.
    local_sup: Point,
}

impl Add<Point> for Extent {
    type Output = Self;

    fn add(self, rhs: Point) -> Extent {
        self.with_minimum(self.minimum + rhs)
    }
}

impl Sub<Point> for Extent {
    type Output = Self;

    fn sub(self, rhs: Point) -> Extent {
        self + (-rhs)
    }
}

impl Extent {
    /// The extent with `minimum` as the least element and `world_sup` as the least upper bound of
    /// the extent, or `minimum + size`.
    pub fn from_min_and_world_supremum(minimum: Point, world_sup: Point) -> Self {
        Self {
            minimum,
            world_sup,
            local_sup: world_sup - minimum,
        }
    }

    /// The extent with `minimum` as the least element and `world_max` as the greatest element.
    pub fn from_min_and_world_max(minimum: Point, world_max: Point) -> Self {
        let world_sup = world_max + [1, 1, 1].into();

        Self {
            minimum,
            world_sup,
            local_sup: world_sup - minimum,
        }
    }

    /// The extent with `minimum` as the least element and `local_max` as the greatest element in
    /// local coordinates (relative to `minimum`).
    pub fn from_minimum_and_local_max(minimum: Point, local_max: Point) -> Self {
        let local_sup = local_max + [1, 1, 1].into();

        Self {
            minimum,
            world_sup: local_sup + minimum,
            local_sup,
        }
    }

    /// The extent with `minimum` as the least element and `local_sup` as the least upper bound
    /// in local coordinates (i.e. the size).
    pub fn from_min_and_local_supremum(minimum: Point, local_sup: Point) -> Self {
        Self {
            minimum,
            world_sup: minimum + local_sup,
            local_sup,
        }
    }

    /// Returns a cube with all dimensions of length `2 * radius + 1`.
    pub fn from_center_and_radius(center: Point, radius: i32) -> Self {
        assert!(radius > 0);

        let minimum = center - [radius; 3].into();
        let local_sup: Point = [2 * radius + 1; 3].into();

        Self::from_min_and_local_supremum(minimum, local_sup)
    }

    /// Returns `self` after growing the size of the extent by `p`.
    pub fn add_to_supremum(&self, p: &Point) -> Self {
        let mut ret = *self;
        ret.local_sup = ret.local_sup + *p;
        ret.world_sup = ret.world_sup + *p;

        ret
    }

    /// Change the size of `self` in-place.
    pub fn set_local_supremum(&mut self, new_sup: Point) {
        *self = Self::from_min_and_local_supremum(self.minimum, new_sup)
    }

    /// Get the center of mass, rounded down to the nearest integer.
    pub fn get_center(&self) -> Point {
        (self.minimum + self.world_sup) / 2
    }

    /// Get the least element.
    pub fn get_minimum(&self) -> Point {
        self.minimum
    }

    /// Get the greatest element.
    pub fn get_world_max(&self) -> Point {
        self.world_sup - [1, 1, 1].into()
    }

    /// Get the greatest element in local coordinates.
    pub fn get_local_max(&self) -> Point {
        self.local_sup - [1, 1, 1].into()
    }

    /// Get the least upper bound in local coordinates (i.e. the size).
    pub fn get_local_supremum(&self) -> &Point {
        &self.local_sup
    }

    /// Get the least upper bound.
    pub fn get_world_supremum(&self) -> &Point {
        &self.world_sup
    }

    /// Translates the entire extent such that `min` is the new minimum.
    pub fn with_minimum(&self, min: Point) -> Self {
        Extent::from_min_and_local_supremum(min, self.local_sup)
    }

    /// Translates the entire extent such that the origin `(0, 0, 0)` is the new minimum.
    pub fn with_minimum_as_origin(&self) -> Self {
        *self - self.minimum
    }

    /// Grows the extent in each of the directions. Positive values grow in the corresponding
    /// direction. Negative values shrink in that direction. Behavior is undefined if the extent
    /// shrinks to be "inside out."
    pub fn directional_grow(&self, grower: &DirectionIndex<i32>) -> Self {
        let positive = grower.positive_point();
        let negative = grower.negative_point();

        let mut ret = *self;
        ret = ret + negative;

        ret.add_to_supremum(&(positive - negative))
    }

    /// Grows (or shrinks for negative `r`) the extent symmetrically in each dimension by `2r`.
    pub fn radial_grow(&self, r: i32) -> Self {
        let grower = DirectionIndex::new([r; 6]);

        self.directional_grow(&grower)
    }

    /// Grows the extent symmetrically in each dimension by `2r`.
    pub fn padded(&self, pad: u32) -> Self {
        self.radial_grow(pad as i32)
    }

    /// Number of lattice points in the extent.
    pub fn volume(&self) -> usize {
        // It's possible for the size of the extent to be negative, but we don't want to return
        // a negative volume.
        (self.local_sup.x.max(0) * self.local_sup.y.max(0) * self.local_sup.z.max(0)) as usize
    }

    /// Returns `true` iff `self` contains no points.
    pub fn is_empty(&self) -> bool {
        self.volume() == 0
    }

    /// Translates `p` from local coordinates to world coordinates.
    pub fn local_point_from_world_point(&self, p: &Point) -> Point {
        *p - self.minimum
    }

    /// Translates `p` from world coordinates to local coordinates.
    pub fn world_point_from_local_point(&self, p: &Point) -> Point {
        *p + self.minimum
    }

    /// Returns true iff `self.get_minimum() + local_point` is an element of `self`.
    pub fn contains_local(&self, local_point: &Point) -> bool {
        Point::new(0, 0, 0) <= *local_point && *local_point < self.local_sup
    }

    /// Returns true iff `world_point` is an element of `self`.
    pub fn contains_world(&self, world_point: &Point) -> bool {
        self.minimum <= *world_point && *world_point < self.world_sup
    }

    /// Returns the overlapping points of `self` and `other`.
    pub fn intersection(&self, other: &Self) -> Self {
        let minimum = self.minimum.join(&other.minimum);
        let world_sup = self.world_sup.meet(&other.world_sup);

        Self::from_min_and_world_supremum(minimum, world_sup)
    }

    /// Returns true `iff` all points in `other` are also in `self`.
    pub fn is_subset(&self, other: &Self) -> bool {
        self.intersection(other) == *self
    }

    /// Get the corner points, i.e. those boundary points which are extreme in each dimension.
    pub fn get_world_corners(&self) -> [Point; 8] {
        let min = self.get_minimum();
        let max = self.get_world_max();

        [
            min,
            [max.x, min.y, min.z].into(),
            [min.x, max.y, min.z].into(),
            [min.x, min.y, max.z].into(),
            [max.x, max.y, min.z].into(),
            [max.x, min.y, max.z].into(),
            [min.x, max.y, max.z].into(),
            max,
        ]
    }

    /// Returns a set of disjoint extents making up the full set of boundary points.
    pub fn get_boundary_extents(&self) -> Vec<Self> {
        let base_area = self.local_sup.x * self.local_sup.y;
        if base_area * self.local_sup.z == 0 {
            return vec![];
        }
        if base_area == 1 {
            return vec![*self];
        }

        let min = self.minimum;
        let max = self.get_local_max();

        let extents = vec![
            Self::from_min_and_world_max((0, 0, 0).into(), (max.x - 1, 0, max.z).into()) + min,
            Self::from_min_and_world_max((0, 1, 0).into(), (0, max.y, max.z).into()) + min,
            Self::from_min_and_world_max((1, max.y, 0).into(), max) + min,
            Self::from_min_and_world_max((max.x, 0, 0).into(), max - (0, 1, 0).into()) + min,
            Self::from_min_and_world_max((1, 1, 0).into(), max - (1, 1, max.z).into()) + min,
            Self::from_min_and_world_max((1, 1, max.z).into(), max - (1, 1, 0).into()) + min,
        ];

        extents.into_iter().filter(|e| !e.is_empty()).collect()
    }

    /// Returns an iterator over the points on the boundary of `self` (those points adjacent to some
    /// point not in `self`).
    pub fn iter_boundary_points(&self) -> impl Iterator<Item = Point> {
        self.get_boundary_extents()
            .into_iter()
            .map(|x| x.into_iter())
            .flatten()
    }

    /// Returns the penetration depths for all normal vectors, such that pushing e1 by the depth
    /// times the corresponding unit normal would always separate the two extents.
    pub fn penetrations(e1: &Self, e2: &Self) -> DirectionIndex<i32> {
        let e1sup = e1.get_world_supremum();
        let e2sup = e2.get_world_supremum();
        let e1min = e1.get_minimum();
        let e2min = e2.get_minimum();

        DirectionIndex {
            values: [
                e1sup.x - e2min.x,
                e2sup.x - e1min.x,
                e1sup.y - e2min.y,
                e2sup.y - e1min.y,
                e1sup.z - e2min.z,
                e2sup.z - e1min.z,
            ],
        }
    }
}

/// Returns the smallest extent containing all of the given points.
pub fn bounding_extent<I>(points: I) -> Extent
where
    I: Iterator<Item = Point>,
{
    let mut min_point = Point::new(std::i32::MAX, std::i32::MAX, std::i32::MAX);
    let mut max_point = Point::new(std::i32::MIN, std::i32::MIN, std::i32::MIN);
    for p in points {
        min_point = min_point.meet(&p);
        max_point = max_point.join(&p);
    }

    Extent::from_min_and_world_max(min_point, max_point)
}

#[derive(Debug)]
pub struct ExtentIterator {
    extent: Extent,
    cursor: Point,
    completed: bool,
}

impl ExtentIterator {
    pub fn new(extent: Extent) -> Self {
        ExtentIterator {
            extent,
            cursor: extent.minimum,
            completed: extent.is_empty(),
        }
    }
}

/// Returns a `Point` for each world coordinate in the extent.
impl Iterator for ExtentIterator {
    type Item = Point;

    fn next(&mut self) -> Option<Point> {
        if self.completed {
            return None;
        }

        let old_cursor = self.cursor;

        self.cursor.z += 1;
        if self.cursor.z == self.extent.world_sup.z {
            self.cursor.z = self.extent.minimum.z;
            self.cursor.y += 1;
            if self.cursor.y == self.extent.world_sup.y {
                self.cursor.y = self.extent.minimum.y;
                self.cursor.x += 1;
                if self.cursor.x == self.extent.world_sup.x {
                    self.completed = true;
                }
            }
        }

        Some(old_cursor)
    }
}

impl IntoIterator for &Extent {
    type Item = Point;
    type IntoIter = ExtentIterator;

    fn into_iter(self) -> Self::IntoIter {
        ExtentIterator::new(*self)
    }
}

impl IntoIterator for Extent {
    type Item = Point;
    type IntoIter = ExtentIterator;

    fn into_iter(self) -> Self::IntoIter {
        ExtentIterator::new(self)
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
    fn test_empty_extent_iterator_generates_nothing() {
        let extent = Extent::from_min_and_world_supremum([0, 0, 0].into(), [1, 0, 2].into());

        let collected: Vec<_> = extent.into_iter().collect();
        assert_eq!(collected, vec![]);
    }

    #[test]
    fn test_nonempty_extent_iterator_generates_correct_points_in_correct_order() {
        let extent = Extent::from_min_and_local_supremum([1, 4, 2].into(), [1, 3, 2].into());

        let collected: Vec<_> = extent.into_iter().collect();
        assert_eq!(
            collected,
            vec![
                [1, 4, 2].into(),
                [1, 4, 3].into(),
                [1, 5, 2].into(),
                [1, 5, 3].into(),
                [1, 6, 2].into(),
                [1, 6, 3].into(),
            ]
        );
    }

    #[test]
    fn test_extent_contains_all_iterator_points() {
        let extent = Extent::from_min_and_world_supremum([0, 0, 0].into(), [1, 3, 2].into());

        for p in &extent {
            assert!(extent.contains_world(&p));
        }
    }

    #[test]
    fn test_extent_does_not_contain_boundary_points() {
        let extent = Extent::from_min_and_local_supremum([1, 1, 1].into(), [1, 3, 2].into());

        // All adjacent to one face of extent.
        let boundaries = vec![
            Extent::from_min_and_world_supremum([0, 0, 0].into(), [1, 10, 10].into()),
            Extent::from_min_and_world_supremum([0, 0, 0].into(), [10, 1, 10].into()),
            Extent::from_min_and_world_supremum([0, 0, 0].into(), [10, 10, 1].into()),
            Extent::from_min_and_local_supremum([2, 1, 1].into(), [10, 10, 10].into()),
            Extent::from_min_and_local_supremum([1, 4, 1].into(), [10, 10, 10].into()),
            Extent::from_min_and_local_supremum([1, 1, 3].into(), [10, 10, 10].into()),
        ];

        for b in &boundaries {
            for p in b {
                assert!(!extent.contains_world(&p));
            }
        }
    }

    #[test]
    fn test_intersection() {
        let e1 = Extent::from_min_and_world_supremum([0, 0, 0].into(), [5, 7, 5].into());
        let e2 = Extent::from_min_and_world_supremum([1, -1, 1].into(), [6, 6, 6].into());

        assert_eq!(
            e1.intersection(&e2),
            Extent::from_min_and_world_supremum([1, 0, 1].into(), [5, 6, 5].into())
        );
    }

    #[test]
    fn test_boundary_of_points() {
        let points: [Point; 5] = [
            (-1, 5, 3).into(),
            (0, 4, 2).into(),
            (8, 0, 1).into(),
            (-4, -8, 0).into(),
            (5, 7, -1).into(),
        ];

        assert_eq!(
            bounding_extent(points.iter().cloned()),
            Extent::from_min_and_world_max((-4, -8, -1).into(), (8, 7, 3).into())
        );
    }

    #[test]
    fn test_get_walls_1x1x1() {
        let e1 = Extent::from_min_and_local_supremum([0, 0, 0].into(), [1, 1, 1].into());

        assert_eq!(e1.get_boundary_extents(), vec![e1]);
    }

    #[test]
    fn test_get_walls_1x1x3() {
        let e1 = Extent::from_min_and_local_supremum([0, 0, 0].into(), [1, 1, 3].into());

        assert_eq!(e1.get_boundary_extents(), vec![e1]);
    }

    #[test]
    fn test_get_walls_2x1x1() {
        let e1 = Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 1, 1].into());

        assert_eq!(
            e1.get_boundary_extents(),
            vec![
                Extent::from_min_and_local_supremum([0, 0, 0].into(), [1, 1, 1].into()),
                Extent::from_min_and_local_supremum([1, 0, 0].into(), [1, 1, 1].into()),
            ]
        );
    }

    #[test]
    fn test_get_walls_2x2x1() {
        let e1 = Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 2, 1].into());

        assert_elements_eq(
            &e1.get_boundary_extents(),
            &vec![
                Extent::from_min_and_local_supremum([0, 0, 0].into(), [1, 1, 1].into()),
                Extent::from_min_and_local_supremum([0, 1, 0].into(), [1, 1, 1].into()),
                Extent::from_min_and_local_supremum([1, 1, 0].into(), [1, 1, 1].into()),
                Extent::from_min_and_local_supremum([1, 0, 0].into(), [1, 1, 1].into()),
            ],
        );
    }

    #[test]
    fn test_get_walls_2x2x2() {
        let e1 = Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 2, 2].into());

        assert_elements_eq(
            &e1.get_boundary_extents(),
            &vec![
                Extent::from_min_and_local_supremum([0, 0, 0].into(), [1, 1, 2].into()),
                Extent::from_min_and_local_supremum([0, 1, 0].into(), [1, 1, 2].into()),
                Extent::from_min_and_local_supremum([1, 1, 0].into(), [1, 1, 2].into()),
                Extent::from_min_and_local_supremum([1, 0, 0].into(), [1, 1, 2].into()),
            ],
        );
    }

    #[test]
    fn test_get_walls_3x3x2() {
        let e1 = Extent::from_min_and_local_supremum([0, 0, 0].into(), [3, 3, 2].into());

        assert_elements_eq(
            &e1.get_boundary_extents(),
            &vec![
                Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 1, 2].into()),
                Extent::from_min_and_local_supremum([0, 1, 0].into(), [1, 2, 2].into()),
                Extent::from_min_and_local_supremum([1, 2, 0].into(), [2, 1, 2].into()),
                Extent::from_min_and_local_supremum([2, 0, 0].into(), [1, 2, 2].into()),
                Extent::from_min_and_local_supremum([1, 1, 0].into(), [1, 1, 1].into()),
                Extent::from_min_and_local_supremum([1, 1, 1].into(), [1, 1, 1].into()),
            ],
        );
    }

    #[test]
    fn test_get_walls_3x3x3() {
        let e1 = Extent::from_min_and_local_supremum([0, 0, 0].into(), [3, 3, 3].into());

        assert_elements_eq(
            &e1.get_boundary_extents(),
            &vec![
                Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 1, 3].into()),
                Extent::from_min_and_local_supremum([0, 1, 0].into(), [1, 2, 3].into()),
                Extent::from_min_and_local_supremum([1, 2, 0].into(), [2, 1, 3].into()),
                Extent::from_min_and_local_supremum([2, 0, 0].into(), [1, 2, 3].into()),
                Extent::from_min_and_local_supremum([1, 1, 0].into(), [1, 1, 1].into()),
                Extent::from_min_and_local_supremum([1, 1, 2].into(), [1, 1, 1].into()),
            ],
        );
    }

    #[test]
    fn test_iter_boundary_points_3x3x3() {
        let e1 = Extent::from_min_and_local_supremum([0, 0, 0].into(), [3, 3, 3].into());

        let expected_points: Vec<_> = e1.into_iter().filter(|p| *p != (1, 1, 1).into()).collect();
        let actual_points: Vec<_> = e1.iter_boundary_points().collect();
        assert_elements_eq(&actual_points, &expected_points);
    }

    #[test]
    fn test_grow_extent_in_negative_directions() {
        let e = Extent::from_min_and_world_max([0, 0, 0].into(), [3, 3, 3].into());

        let grown_e = e.directional_grow(&DirectionIndex {
            values: [1, 0, 1, 0, 1, 0],
        });

        assert_eq!(
            grown_e,
            Extent::from_min_and_world_max([-1, -1, -1].into(), [3, 3, 3].into(),)
        );
    }

    #[test]
    fn test_grow_extent_in_positive_directions() {
        let e = Extent::from_min_and_world_max([0, 0, 0].into(), [3, 3, 3].into());

        let grown_e = e.directional_grow(&DirectionIndex {
            values: [0, 1, 0, 1, 0, 1],
        });

        assert_eq!(
            grown_e,
            Extent::from_min_and_world_max([0, 0, 0].into(), [4, 4, 4].into(),)
        );
    }

    #[test]
    fn test_shrink_extent_in_negative_directions() {
        let e = Extent::from_min_and_world_max([0, 0, 0].into(), [3, 3, 3].into());

        let shrunk_e = e.directional_grow(&DirectionIndex {
            values: [-1, 0, -1, 0, -1, 0],
        });

        assert_eq!(
            shrunk_e,
            Extent::from_min_and_world_max([1, 1, 1].into(), [3, 3, 3].into(),)
        );
    }

    #[test]
    fn test_shrink_extent_in_positive_directions() {
        let e = Extent::from_min_and_world_max([0, 0, 0].into(), [3, 3, 3].into());

        let shrunk_e = e.directional_grow(&DirectionIndex {
            values: [0, -1, 0, -1, 0, -1],
        });

        assert_eq!(
            shrunk_e,
            Extent::from_min_and_world_max([0, 0, 0].into(), [2, 2, 2].into(),)
        );
    }
}
