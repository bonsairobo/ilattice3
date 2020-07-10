use crate::{copy_extent, Extent, GetWorld, Point, VecLatticeMap};

pub struct FnLatticeMap<F> {
    f: F,
}

impl<F> FnLatticeMap<F> {
    pub fn new(f: F) -> Self {
        FnLatticeMap { f }
    }

    /// Samples `self` at all points in `extent`, and stores the values into a `VecLatticeMap`.
    pub fn render<T>(&self, extent: &Extent) -> VecLatticeMap<T>
    where
        T: Clone + Default,
        F: Fn(&Point) -> T,
    {
        let mut new_map = VecLatticeMap::fill(*extent, T::default());
        copy_extent(self, &mut new_map, extent);

        new_map
    }
}

impl<T, F> GetWorld<T> for FnLatticeMap<F>
where
    F: Fn(&Point) -> T,
{
    fn get_world(&self, p: &Point) -> T {
        (self.f)(p)
    }
}
