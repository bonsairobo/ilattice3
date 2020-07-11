use crate::{GetWorld, Point};

pub struct FnLatticeMap<F> {
    f: F,
}

impl<F> FnLatticeMap<F> {
    pub fn new(f: F) -> Self {
        FnLatticeMap { f }
    }
}

impl<T, F> GetWorld for FnLatticeMap<F>
where
    F: Fn(&Point) -> T,
{
    type Data = T;
    fn get_world(&self, p: &Point) -> T {
        (self.f)(p)
    }
}
