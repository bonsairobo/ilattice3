//! A bunch of common traits for the various kinds of lattice maps.

use crate::{extent::ExtentIterator, GetExtent, HasIndexer, Indexer, Point};

use std::borrow::Borrow;

/// Returns the data at index `i`.
pub trait GetLinear {
    type Data;
    fn get_linear(&self, i: usize) -> Self::Data;
}

/// Returns a reference to the data at index `i`.
pub trait GetLinearRef {
    type Data;
    fn get_linear_ref(&self, i: usize) -> &Self::Data;
}

/// Returns a mutable reference to the data at index `i`.
pub trait GetLinearRefMut {
    type Data;
    fn get_linear_ref_mut(&mut self, i: usize) -> &mut Self::Data;
}

/// Returns the data at point `p` in lattice-local coordinates, i.e. relative to some
/// fixed origin.
pub trait GetLocal {
    type Data;
    fn get_local(&self, p: &Point) -> Self::Data;
}

/// Returns a reference to the data at point `p` in lattice-local coordinates, i.e. relative to some
/// fixed origin.
pub trait GetLocalRef {
    type Data;
    fn get_local_ref(&self, p: &Point) -> &Self::Data;
}

/// Returns a mutable reference to the data at point `p` in lattice-local coordinates, i.e. relative
/// to some fixed origin.
pub trait GetLocalRefMut {
    type Data;
    fn get_local_ref_mut(&mut self, p: &Point) -> &mut Self::Data;
}

/// Returns the data at point `p` in world coordinates, i.e. the origin is [0, 0, 0].
pub trait GetWorld {
    type Data;
    fn get_world(&self, p: &Point) -> Self::Data;
}

/// Returns a reference to the data at point `p` in world coordinates, i.e. the origin is [0, 0, 0].
pub trait GetWorldRef {
    type Data;
    fn get_world_ref(&self, p: &Point) -> &Self::Data;
}

/// Returns a mutable reference to the data at point `p` in world coordinates, i.e. the origin is
/// [0, 0, 0].
pub trait GetWorldRefMut {
    type Data;
    fn get_world_ref_mut(&mut self, p: &Point) -> &mut Self::Data;
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

pub trait MaybeGetWorld {
    type Data;
    fn maybe_get_world(&self, p: &Point) -> Option<Self::Data>;
}

pub trait MaybeGetWorldRef {
    type Data;
    fn maybe_get_world_ref(&self, p: &Point) -> Option<&Self::Data>;
}

pub trait MaybeGetWorldRefMut {
    type Data;
    fn maybe_get_world_ref_mut(&mut self, p: &Point) -> Option<&mut Self::Data>;
}

//   _     _             _        _
//  | |__ | | __ _ _ __ | | _____| |_ ___
//  | '_ \| |/ _` | '_ \| |/ / _ \ __/ __|
//  | |_) | | (_| | | | |   <  __/ |_\__ \
//  |_.__/|_|\__,_|_| |_|_|\_\___|\__|___/

impl<M, T, I> GetLocal for M
where
    M: GetExtent + GetLinear<Data = T> + HasIndexer<Indexer = I>,
    I: Indexer,
{
    type Data = T;

    fn get_local(&self, p: &Point) -> Self::Data {
        let s = self.get_extent().get_local_supremum();

        self.get_linear(I::index_from_local_point(s, p))
    }
}

impl<M, T, I> GetLocalRef for M
where
    M: GetExtent + GetLinearRef<Data = T> + HasIndexer<Indexer = I>,
    I: Indexer,
{
    type Data = T;

    fn get_local_ref(&self, p: &Point) -> &Self::Data {
        let s = *self.get_extent().get_local_supremum();

        self.get_linear_ref(I::index_from_local_point(&s, p))
    }
}

impl<M, T, I> GetLocalRefMut for M
where
    M: GetExtent + GetLinearRefMut<Data = T> + HasIndexer<Indexer = I>,
    I: Indexer,
{
    type Data = T;

    fn get_local_ref_mut(&mut self, p: &Point) -> &mut Self::Data {
        let s = *self.get_extent().get_local_supremum();

        self.get_linear_ref_mut(I::index_from_local_point(&s, p))
    }
}

impl<M, T> GetWorld for M
where
    M: GetLocal<Data = T> + GetExtent,
{
    type Data = T;

    fn get_world(&self, p: &Point) -> Self::Data {
        let p_local = self.get_extent().local_point_from_world_point(p);

        self.get_local(&p_local)
    }
}

impl<M, T> GetWorldRef for M
where
    M: GetLocalRef<Data = T> + GetExtent,
{
    type Data = T;

    fn get_world_ref(&self, p: &Point) -> &T {
        let p_local = self.get_extent().local_point_from_world_point(p);

        self.get_local_ref(&p_local)
    }
}

impl<M, T> GetWorldRefMut for M
where
    M: GetLocalRefMut<Data = T> + GetExtent,
{
    type Data = T;

    fn get_world_ref_mut(&mut self, p: &Point) -> &mut T {
        let p_local = self.get_extent().local_point_from_world_point(p);

        self.get_local_ref_mut(&p_local)
    }
}

impl<'a, M, T> GetWorldBorrowable<'a, T, T> for M
where
    M: GetWorld<Data = T>,
{
    fn get_world_borrowable<'b: 'a>(&'b self, p: &Point) -> T {
        self.get_world(p)
    }
}

impl<'a, M, T> GetWorldBorrowable<'a, T, &'a T> for M
where
    M: GetWorldRef<Data = T>,
{
    fn get_world_borrowable<'b: 'a>(&'b self, p: &Point) -> &'a T {
        self.get_world_ref(p)
    }
}

impl<M, T> MaybeGetWorld for M
where
    M: GetWorld<Data = T> + GetExtent,
{
    type Data = T;

    fn maybe_get_world(&self, p: &Point) -> Option<T> {
        if self.get_extent().contains_world(p) {
            Some(self.get_world(p))
        } else {
            None
        }
    }
}

impl<M, T> MaybeGetWorldRef for M
where
    M: GetWorldRef<Data = T> + GetExtent,
{
    type Data = T;

    fn maybe_get_world_ref(&self, p: &Point) -> Option<&T> {
        if self.get_extent().contains_world(p) {
            Some(self.get_world_ref(p))
        } else {
            None
        }
    }
}

impl<M, T> MaybeGetWorldRefMut for M
where
    M: GetWorldRefMut<Data = T> + GetExtent,
{
    type Data = T;

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
    V: GetWorldRef<Data = T>,
    T: 'a,
{
    type Item = (Point, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        self.extent_iter
            .next()
            .map(|p| (p, self.map.get_world_ref(&p)))
    }
}

/// Wraps any `V: GetExtent + GetWorldRef<T>` to make it `IntoIterator`.
pub struct LatticeMapIter<V>(pub V);

// This impl requires that `GetWorldRef` has an associated `Data` type, since otherwise `T` would be
// unconstrained.
impl<'a, V, T> IntoIterator for LatticeMapIter<&'a V>
where
    V: GetExtent + GetWorldRef<Data = T>,
    T: 'a,
{
    type Item = (Point, &'a T);
    type IntoIter = LatticeMapKeyValIterator<'a, V, T>;

    fn into_iter(self) -> Self::IntoIter {
        LatticeMapKeyValIterator::new(&self.0, self.0.get_extent().into_iter())
    }
}
