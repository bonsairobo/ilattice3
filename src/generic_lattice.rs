use crate::{Extent, Point};

/// Get the `Extent` for some `Lattice`-like type.
pub trait GetExtent {
    fn get_extent(&self) -> Extent;
}

/// Returns a reference to the data at point `p` in lattice-local coordinates, i.e. the minimum of
/// the extent is [0, 0, 0].
pub trait GetLocal<T> {
    fn get_local(&self, p: &Point) -> &T;
}

/// Returns a mutable reference to the data at point `p` in lattice-local coordinates, i.e. the
/// minimum of the extent is [0, 0, 0].
pub trait GetLocalMut<T> {
    fn get_mut_local(&mut self, p: &Point) -> &mut T;
}

/// Returns a reference to the data at point `p` in world coordinates.
pub trait GetWorld<T> {
    fn get_world(&self, p: &Point) -> &T;

    fn all<F>(&self, extent: &Extent, f: F) -> bool
    where
        F: Fn(&T) -> bool,
    {
        for p in extent {
            if !f(self.get_world(&p)) {
                return false;
            }
        }

        true
    }

    fn some<F>(&self, extent: &Extent, f: F) -> bool
    where
        F: Fn(&T) -> bool,
    {
        for p in extent {
            if f(self.get_world(&p)) {
                return true;
            }
        }

        false
    }
}

/// Returns a mutable reference to the data at point `p` in world coordinates.
pub trait GetWorldMut<T> {
    fn get_mut_world(&mut self, p: &Point) -> &mut T;
}

pub trait MaybeGetWorld<T> {
    fn maybe_get_world(&self, p: &Point) -> Option<&T>;
}

pub trait MaybeGetWorldMut<T> {
    fn maybe_get_world_mut(&mut self, p: &Point) -> Option<&mut T>;
}

impl<L, T> GetWorld<T> for L
where
    L: GetLocal<T> + GetExtent,
{
    fn get_world(&self, p: &Point) -> &T {
        let p_local = self.get_extent().local_point_from_world_point(p);

        self.get_local(&p_local)
    }
}

impl<L, T> GetWorldMut<T> for L
where
    L: GetLocalMut<T> + GetExtent,
{
    fn get_mut_world(&mut self, p: &Point) -> &mut T {
        let p_local = self.get_extent().local_point_from_world_point(p);

        self.get_mut_local(&p_local)
    }
}

pub fn map_extent<T, R, S, D, F>(src: &S, dst: &mut D, extent: &Extent, f: F)
where
    F: Fn(&T) -> R,
    S: GetWorld<T>,
    D: GetWorldMut<R>,
{
    for p in extent {
        *dst.get_mut_world(&p) = f(src.get_world(&p));
    }
}

pub fn fill_extent<T: Clone, D: GetWorldMut<T>>(dst: &mut D, extent: &Extent, val: T) {
    for p in extent {
        *dst.get_mut_world(&p) = val.clone();
    }
}

pub fn copy_extent_to_position<T, R, S, D>(
    src: &S,
    dst: &mut D,
    dst_position: &Point,
    extent: &Extent,
) where
    T: Clone,
    R: From<T>,
    S: GetWorld<T>,
    D: GetWorldMut<R>,
{
    for p in extent {
        let p_dst = *dst_position + p - extent.get_minimum();
        *dst.get_mut_world(&p_dst) = src.get_world(&p).clone().into();
    }
}

pub fn copy_extent<T, R, S, D>(src: &S, dst: &mut D, extent: &Extent)
where
    T: Clone,
    R: From<T>,
    S: GetWorld<T>,
    D: GetWorldMut<R>,
{
    copy_extent_to_position(src, dst, &extent.get_minimum(), extent)
}
