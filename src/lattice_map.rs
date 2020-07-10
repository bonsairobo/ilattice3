use crate::{extent::ExtentIterator, GetExtent, Point};

use std::borrow::Borrow;

/// Returns the data at point `p` in lattice-local coordinates, i.e. relative to some
/// fixed origin.
pub trait GetLocal<T> {
    fn get_local(&self, p: &Point) -> T;
}

/// Returns a reference to the data at point `p` in lattice-local coordinates, i.e. relative to some
/// fixed origin.
pub trait GetLocalRef<T> {
    fn get_local_ref(&self, p: &Point) -> &T;
}

/// Returns a mutable reference to the data at point `p` in lattice-local coordinates, i.e. relative
/// to some fixed origin.
pub trait GetLocalRefMut<T> {
    fn get_local_ref_mut(&mut self, p: &Point) -> &mut T;
}

/// Returns the data at point `p` in world coordinates, i.e. the origin is [0, 0, 0].
pub trait GetWorld<T> {
    fn get_world(&self, p: &Point) -> T;
}

/// Returns a reference to the data at point `p` in world coordinates, i.e. the origin is [0, 0, 0].
pub trait GetWorldRef<T> {
    fn get_world_ref(&self, p: &Point) -> &T;
}

/// Returns a mutable reference to the data at point `p` in world coordinates, i.e. the origin is
/// [0, 0, 0].
pub trait GetWorldRefMut<T> {
    fn get_world_ref_mut(&mut self, p: &Point) -> &mut T;
}

/// Returns some type that can be borrowed at the point `p`. Mostly an implementation detail for
/// functions that need to be generic over lattice maps that implement `GetWorld` or `GetWorldRef`.
pub trait GetWorldBorrowable<'a, T, B>
where
    B: Borrow<T>,
{
    fn get_world_borrowable<'b>(&'b self, p: &Point) -> B
    where
        'b: 'a;
}

pub trait MaybeGetWorld<T> {
    fn maybe_get_world(&self, p: &Point) -> Option<T>;
}

pub trait MaybeGetWorldRef<T> {
    fn maybe_get_world_ref(&self, p: &Point) -> Option<&T>;
}

pub trait MaybeGetWorldRefMut<T> {
    fn maybe_get_world_ref_mut(&mut self, p: &Point) -> Option<&mut T>;
}

//   _     _             _        _
//  | |__ | | __ _ _ __ | | _____| |_ ___
//  | '_ \| |/ _` | '_ \| |/ / _ \ __/ __|
//  | |_) | | (_| | | | |   <  __/ |_\__ \
//  |_.__/|_|\__,_|_| |_|_|\_\___|\__|___/

impl<M, T> GetWorld<T> for M
where
    M: GetLocal<T> + GetExtent,
{
    fn get_world(&self, p: &Point) -> T {
        let p_local = self.get_extent().local_point_from_world_point(p);

        self.get_local(&p_local)
    }
}

impl<M, T> GetWorldRef<T> for M
where
    M: GetLocalRef<T> + GetExtent,
{
    fn get_world_ref(&self, p: &Point) -> &T {
        let p_local = self.get_extent().local_point_from_world_point(p);

        self.get_local_ref(&p_local)
    }
}

impl<M, T> GetWorldRefMut<T> for M
where
    M: GetLocalRefMut<T> + GetExtent,
{
    fn get_world_ref_mut(&mut self, p: &Point) -> &mut T {
        let p_local = self.get_extent().local_point_from_world_point(p);

        self.get_local_ref_mut(&p_local)
    }
}

impl<'a, M, T> GetWorldBorrowable<'a, T, T> for M
where
    M: GetWorld<T>,
{
    fn get_world_borrowable<'b: 'a>(&'b self, p: &Point) -> T {
        self.get_world(p)
    }
}

impl<'a, M, T> GetWorldBorrowable<'a, T, &'a T> for M
where
    M: GetWorldRef<T>,
{
    fn get_world_borrowable<'b: 'a>(&'b self, p: &Point) -> &'a T {
        self.get_world_ref(p)
    }
}

impl<M, T> MaybeGetWorld<T> for M
where
    M: GetWorld<T> + GetExtent,
{
    fn maybe_get_world(&self, p: &Point) -> Option<T> {
        if self.get_extent().contains_world(p) {
            Some(self.get_world(p))
        } else {
            None
        }
    }
}

impl<M, T> MaybeGetWorldRef<T> for M
where
    M: GetWorldRef<T> + GetExtent,
{
    fn maybe_get_world_ref(&self, p: &Point) -> Option<&T> {
        if self.get_extent().contains_world(p) {
            Some(self.get_world_ref(p))
        } else {
            None
        }
    }
}

impl<M, T> MaybeGetWorldRefMut<T> for M
where
    M: GetWorldRefMut<T> + GetExtent,
{
    fn maybe_get_world_ref_mut(&mut self, p: &Point) -> Option<&mut T> {
        if self.get_extent().contains_world(p) {
            Some(self.get_world_ref_mut(p))
        } else {
            None
        }
    }
}

// There are more blanket impls I would like once Rust has specialization:
// - MaybeGetX for any GetX that just always returns Some
// - GetWorld<T> for any GetWorldRef<T> where T: Clone

/// An iterator over the points and data in a map.
pub struct LatticeMapKeyValIterator<'a, V, T> {
    map: &'a V,
    extent_iter: ExtentIterator,
    marker: std::marker::PhantomData<T>,
}

impl<'a, V, T> LatticeMapKeyValIterator<'a, V, T> {
    pub fn new(map: &'a V, extent_iter: ExtentIterator) -> Self {
        Self {
            map,
            extent_iter,
            marker: Default::default(),
        }
    }
}

impl<'a, V, T> Iterator for LatticeMapKeyValIterator<'a, V, T>
where
    V: GetWorldRef<T>,
    T: 'a,
{
    type Item = (Point, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        self.extent_iter
            .next()
            .map(|p| (p, self.map.get_world_ref(&p)))
    }
}
