use crate::Point;

pub trait PeriodicIndexer: Indexer {}

pub trait Indexer: Clone {
    /// `s` is the local strict supremum of an extent. `p` is a local point.
    fn index_from_local_point(s: &Point, p: &Point) -> usize;

    /// `s` is the local strict supremum of an extent. `index` is a linear index.
    fn local_point_from_index(s: &Point, index: usize) -> Point;
}

/// Most `Indexer`s should not require state to be instantiated.
pub trait StatelessIndexer: Indexer {
    fn new() -> Self;
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct YLevelsIndexer;

impl StatelessIndexer for YLevelsIndexer {
    fn new() -> Self {
        YLevelsIndexer {}
    }
}

impl Indexer for YLevelsIndexer {
    fn index_from_local_point(s: &Point, p: &Point) -> usize {
        // This scheme is chosen for ease of hand-crafting voxel maps.
        (p.y * s.x * s.z + p.z * s.x + p.x) as usize
    }

    fn local_point_from_index(s: &Point, index: usize) -> Point {
        assert!(index <= std::i32::MAX as usize);
        let index = index as i32;
        let xz_area = s.x * s.z;
        let y = index / xz_area;
        let rem = index - y * xz_area;
        let z = rem / s.x;
        let rem = rem - z * s.x;
        let x = rem;

        Point::new(x, y, z)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PeriodicYLevelsIndexer;

impl PeriodicIndexer for PeriodicYLevelsIndexer {}

impl StatelessIndexer for PeriodicYLevelsIndexer {
    fn new() -> Self {
        PeriodicYLevelsIndexer {}
    }
}

impl Indexer for PeriodicYLevelsIndexer {
    fn index_from_local_point(s: &Point, p: &Point) -> usize {
        let px = p.x.rem_euclid(s.x);
        let py = p.y.rem_euclid(s.y);
        let pz = p.z.rem_euclid(s.z);

        (py * s.x * s.z + pz * s.x + px) as usize
    }

    /// Returns the canonical point which is contained in the extent.
    fn local_point_from_index(s: &Point, index: usize) -> Point {
        YLevelsIndexer::local_point_from_index(s, index)
    }
}