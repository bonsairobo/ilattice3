use crate::{copy_extent, prelude::*, Extent, Transform, YLevelsIndexer};

use compressible_map::{Compressible, Decompressible};
use serde::{de::DeserializeOwned, Deserialize, Serialize};
use std::marker::PhantomData;

/// A map from points in an extent to some kind of data `T`, stored as a `Vec<T>`.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct VecLatticeMap<T, I = YLevelsIndexer> {
    // Works if I: Default.
    #[serde(skip_deserializing)]
    #[serde(skip_serializing)]
    indexer: PhantomData<I>,

    extent: Extent,
    values: Vec<T>,
}

impl<T, I> GetExtent for VecLatticeMap<T, I> {
    fn get_extent(&self) -> &Extent {
        &self.extent
    }
}

impl<T, I> HasIndexer for VecLatticeMap<T, I>
where
    I: Indexer,
{
    type Indexer = I;
}

impl<T, I: Indexer> GetLinear for VecLatticeMap<T, I>
where
    T: Clone,
{
    type Data = T;

    fn get_linear(&self, i: usize) -> T {
        self.values[i].clone()
    }
}

impl<T, I: Indexer> GetLinearRef for VecLatticeMap<T, I> {
    type Data = T;

    fn get_linear_ref(&self, i: usize) -> &T {
        &self.values[i]
    }
}

impl<T, I: Indexer> GetLinearRefMut for VecLatticeMap<T, I> {
    type Data = T;

    fn get_linear_ref_mut(&mut self, i: usize) -> &mut T {
        &mut self.values[i]
    }
}

impl<T, I: Indexer> VecLatticeMap<T, I> {
    /// Creates a new `VecLatticeMap` with the given `values`, which are assumed to be ordered
    /// according to a newly created `indexer: I`.
    pub fn new(extent: Extent, values: Vec<T>) -> Self {
        VecLatticeMap {
            extent,
            values,
            indexer: Default::default(),
        }
    }

    /// Same as `Self::new`, but with minimum at the origin, and dimensions (or supremum) `sup`.
    pub fn new_at_origin(sup: Point, values: Vec<T>) -> Self {
        let extent = Extent::from_min_and_world_supremum([0, 0, 0].into(), sup);

        Self::new(extent, values)
    }

    /// Get a linear index from the point `p` in local coordinates using the indexer of `self`.
    pub fn index_from_local_point(&self, p: &Point) -> usize {
        let local_sup = self.extent.get_local_supremum();

        I::index_from_local_point(&local_sup, p)
    }

    /// Get a point `p` from a linear `index` using the indexer of `self`.
    pub fn local_point_from_index(&self, index: usize) -> Point {
        let local_sup = self.extent.get_local_supremum();

        I::local_point_from_index(&local_sup, index)
    }

    /// Get a linear index from the point `p` in world coordinates using the indexer of `self`.
    pub fn index_from_world_point(&self, p: &Point) -> usize {
        self.index_from_local_point(&self.extent.local_point_from_world_point(p))
    }

    /// Translate the entire extent by `delta`.
    pub fn translate(&mut self, delta: &Point) {
        self.extent = self.extent + *delta;
    }

    /// Move the minimum of the extent to `new_min`.
    pub fn set_minimum(&mut self, new_min: &Point) {
        self.extent = self.extent.with_minimum(*new_min);
    }

    /// Create a new lattice map by applying `f` pointwise.
    pub fn map<F, S>(&self, f: F) -> VecLatticeMap<S, I>
    where
        F: Fn(&T) -> S,
    {
        VecLatticeMap::new(*self.get_extent(), self.values.iter().map(f).collect())
    }
}

impl<T: Clone, I: Indexer> VecLatticeMap<T, I> {
    /// Map every point by `tfm`. This function will assert `tfm.is_octahedral` in debug mode.
    pub fn apply_octahedral_transform(&self, tfm: &Transform) -> Self {
        debug_assert!(tfm.is_octahedral());

        let extent = self.get_extent();
        let volume = extent.volume();

        let tfm_extent = tfm.apply_to_extent(&extent);

        let mut new_values = Vec::with_capacity(volume);
        unsafe {
            new_values.set_len(volume);
        }
        let mut tfm_map = Self::new(tfm_extent, new_values);

        // PERF: this is not the most efficient, but it is very simple.
        for p in extent {
            let tfm_p = tfm.apply_to_point(&p);
            *tfm_map.get_world_ref_mut(&tfm_p) = self.get_world(&p);
        }

        tfm_map
    }

    /// Returns a vec of the data in `extent`, ordered linearly by `I: Indexer`. A map
    /// can be recreated from the vec using `VecLatticeMap::<T, I>::new_with_indexer`.
    pub fn linearize_extent(&self, extent: &Extent) -> Vec<T> {
        let num_elements = extent.volume();
        let mut data = Vec::with_capacity(num_elements);
        unsafe {
            data.set_len(num_elements);
        }
        for p in extent {
            let i = I::index_from_local_point(
                extent.get_local_supremum(),
                &extent.local_point_from_world_point(&p),
            );
            data[i] = self.get_world(&p);
        }

        data
    }

    /// Copy all values in `extent` to a new `VecLatticeMap` of the same extent.
    pub fn copy_extent_into_new_map(&self, extent: &Extent) -> Self {
        let volume = extent.volume();
        let mut values = Vec::with_capacity(volume);
        unsafe {
            values.set_len(volume);
        }
        let mut copy = VecLatticeMap::new(*extent, values);
        copy_extent(self, &mut copy, extent);

        copy
    }

    /// Set the value at every point in `extent` to `init_val`.
    pub fn fill(extent: Extent, init_val: T) -> Self {
        VecLatticeMap {
            extent,
            values: vec![init_val; extent.volume()],
            indexer: Default::default(),
        }
    }

    /// Creates a new `VecLatticeMap` by copying the values from another `map` at all points in
    /// `extent`.
    pub fn copy_from_map<V>(map: &V, extent: &Extent) -> Self
    where
        V: GetWorld<Data = T>,
    {
        let volume = extent.volume();
        let mut values = Vec::with_capacity(volume);
        unsafe {
            values.set_len(volume);
        }
        let mut copy = VecLatticeMap::new(*extent, values);
        copy_extent(map, &mut copy, extent);

        copy
    }
}

#[derive(Clone)]
pub struct PortableLZ4 {
    pub level: u32,
}

/// A compressed `VecLatticeMap` that can be stored and decompressed safely on any platform.
#[derive(Clone)]
pub struct PortableCompressedVecLatticeMap<T, I> {
    compressed_bytes: Vec<u8>,
    marker: std::marker::PhantomData<(T, I)>,
}

impl<T, I> Decompressible<PortableLZ4> for PortableCompressedVecLatticeMap<T, I>
where
    T: DeserializeOwned + Serialize,
{
    type Decompressed = VecLatticeMap<T, I>;

    fn decompress(&self) -> Self::Decompressed {
        let mut decoder = lz4::Decoder::new(self.compressed_bytes.as_slice()).unwrap();
        let mut decompressed_bytes = Vec::new();
        std::io::copy(&mut decoder, &mut decompressed_bytes).unwrap();

        bincode::deserialize(decompressed_bytes.as_slice()).unwrap()
    }
}

impl<T, I> Compressible<PortableLZ4> for VecLatticeMap<T, I>
where
    T: DeserializeOwned + Serialize,
{
    type Compressed = PortableCompressedVecLatticeMap<T, I>;

    /// Compress the map in-memory using the LZ4 algorithm. The result should be safe to decompress
    /// on another platform.
    ///
    /// TODO: maybe this should use a different algorithm with better compression ratio, if it's
    /// intended for long-term storage
    fn compress(&self, params: PortableLZ4) -> Self::Compressed {
        let serialized_bytes = bincode::serialize(&self).unwrap();

        let mut compressed_bytes = Vec::new();
        let mut encoder = lz4::EncoderBuilder::new()
            .level(params.level)
            .build(&mut compressed_bytes)
            .unwrap();

        std::io::copy(&mut std::io::Cursor::new(serialized_bytes), &mut encoder).unwrap();
        let (_output, _result) = encoder.finish();

        PortableCompressedVecLatticeMap {
            compressed_bytes,
            marker: Default::default(),
        }
    }
}

#[derive(Clone)]
pub struct FastLZ4 {
    pub level: u32,
}

/// A compressed `VecLatticeMap` that decompresses quickly, but only on the same platform where it
/// was compressed.
#[derive(Clone)]
pub struct FastCompressedVecLatticeMap<T, I> {
    compressed_bytes: Vec<u8>,
    extent: Extent,
    marker: std::marker::PhantomData<(T, I)>,
}

impl<T, I> FastCompressedVecLatticeMap<T, I> {
    pub fn get_extent(&self) -> &Extent {
        &self.extent
    }
}

impl<T, I> Decompressible<FastLZ4> for FastCompressedVecLatticeMap<T, I>
where
    I: Indexer,
{
    type Decompressed = VecLatticeMap<T, I>;

    fn decompress(&self) -> Self::Decompressed {
        let volume = self.extent.volume();

        let mut decoder = lz4::Decoder::new(self.compressed_bytes.as_slice()).unwrap();
        let mut decompressed_values: Vec<T> = Vec::with_capacity(volume);
        unsafe { decompressed_values.set_len(volume) };
        let mut decompressed_slice = unsafe {
            std::slice::from_raw_parts_mut(
                decompressed_values.as_mut_ptr() as *mut u8,
                volume * std::mem::size_of::<T>(),
            )
        };
        std::io::copy(&mut decoder, &mut decompressed_slice).unwrap();

        VecLatticeMap::new(self.extent, decompressed_values)
    }
}

impl<T, I> Compressible<FastLZ4> for VecLatticeMap<T, I>
where
    I: Indexer,
{
    type Compressed = FastCompressedVecLatticeMap<T, I>;

    /// Compress the map in-memory using the LZ4 algorithm.
    ///
    /// WARNING: For performance, this reinterprets the inner vector as a byte slice without
    /// accounting for endianness. This is not compatible across platforms.
    fn compress(&self, params: FastLZ4) -> FastCompressedVecLatticeMap<T, I> {
        let mut compressed_bytes = Vec::new();
        let values_slice: &[u8] = unsafe {
            std::slice::from_raw_parts(
                self.values.as_ptr() as *const u8,
                self.values.len() * std::mem::size_of::<T>(),
            )
        };
        let mut encoder = lz4::EncoderBuilder::new()
            .level(params.level)
            .build(&mut compressed_bytes)
            .unwrap();

        std::io::copy(&mut std::io::Cursor::new(values_slice), &mut encoder).unwrap();
        let (_output, _result) = encoder.finish();

        FastCompressedVecLatticeMap {
            extent: self.extent,
            compressed_bytes,
            marker: Default::default(),
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
    use crate::PeriodicYLevelsIndexer;

    #[test]
    fn test_periodic_indexer() {
        let extent = Extent::from_min_and_local_supremum([-1, -1, -1].into(), [3, 3, 3].into());
        let mut map = VecLatticeMap::<_, PeriodicYLevelsIndexer>::fill(extent, (0, 0, 0));
        for p in &extent {
            *map.get_world_ref_mut(&p) = (p.x, p.y, p.z);
        }

        assert_eq!(map.get_world(&[-1, -1, -1].into()), (-1, -1, -1));

        assert_eq!(map.get_world(&[-2, -1, -1].into()), (1, -1, -1));
        assert_eq!(map.get_world(&[-1, -2, -1].into()), (-1, 1, -1));
        assert_eq!(map.get_world(&[-1, -1, -2].into()), (-1, -1, 1));

        assert_eq!(map.get_world(&[-3, -1, -1].into()), (0, -1, -1));
        assert_eq!(map.get_world(&[-1, -3, -1].into()), (-1, 0, -1));
        assert_eq!(map.get_world(&[-1, -1, -3].into()), (-1, -1, 0));
    }

    #[test]
    fn test_local_point_from_index() {
        let sup = Point::new(10, 20, 30);
        let test_points = [
            Point::new(0, 0, 0),
            Point::new(1, 1, 1),
            Point::new(1, 2, 3),
            Point::new(4, 3, 2),
            Point::new(9, 9, 10),
            Point::new(9, 19, 29),
        ];

        for p in test_points.iter() {
            assert_eq!(
                *p,
                YLevelsIndexer::local_point_from_index(
                    &sup,
                    YLevelsIndexer::index_from_local_point(&sup, p)
                ),
            );
        }
    }

    #[test]
    fn test_serialized_extent_back_to_map() {
        let extent = Extent::from_min_and_local_supremum([-1, -1, -1].into(), [3, 3, 3].into());
        let mut map = VecLatticeMap::<_, YLevelsIndexer>::fill(extent, (0, 0, 0));
        for p in &extent {
            *map.get_world_ref_mut(&p) = (p.x, p.y, p.z);
        }
        let orig_map = map.clone();

        let serial_extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 2, 2].into());
        let serialized = map.linearize_extent(&serial_extent);
        let serial_map = VecLatticeMap::<_, YLevelsIndexer>::new(serial_extent, serialized);
        copy_extent(&serial_map, &mut map, &serial_extent);

        assert_eq!(orig_map, map);
    }

    #[test]
    fn test_apply_identity() {
        let matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let tfm = Transform { matrix };

        let extent = Extent::from_min_and_local_supremum([-1, -1, -1].into(), [3, 3, 3].into());
        let mut map = VecLatticeMap::<_, YLevelsIndexer>::fill(extent, (0, 0, 0));
        for p in &extent {
            *map.get_world_ref_mut(&p) = (p.x, p.y, p.z);
        }

        let orig_map = map.clone();
        let tfm_map = map.apply_octahedral_transform(&tfm);

        assert_eq!(orig_map, tfm_map);
    }

    #[test]
    fn test_apply_x_reflection() {
        let matrix = [[-1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let tfm = Transform { matrix };

        let orig_extent =
            Extent::from_min_and_local_supremum([-2, -2, -2].into(), [3, 3, 3].into());
        let mut orig_map = VecLatticeMap::<_, YLevelsIndexer>::fill(orig_extent, (0, 0, 0));
        for p in &orig_extent {
            *orig_map.get_world_ref_mut(&p) = (p.x, p.y, p.z);
        }

        let tfm_map = orig_map.apply_octahedral_transform(&tfm);

        let manual_tfm_extent =
            Extent::from_min_and_local_supremum([0, -2, -2].into(), [3, 3, 3].into());
        let mut manual_tfm_map = VecLatticeMap::fill(manual_tfm_extent, (0, 0, 0));
        for p in &manual_tfm_extent {
            *manual_tfm_map.get_world_ref_mut(&p) = (-p.x, p.y, p.z);
        }

        assert_eq!(tfm_map, manual_tfm_map);
    }

    #[test]
    fn test_apply_rotation() {
        let matrix = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]];
        let tfm = Transform { matrix };

        let orig_extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), [1, 2, 3].into());
        let mut orig_map = VecLatticeMap::<_, YLevelsIndexer>::fill(orig_extent, (0, 0, 0));
        for p in &orig_extent {
            *orig_map.get_world_ref_mut(&p) = (p.x, p.y, p.z);
        }

        let tfm_map = orig_map.apply_octahedral_transform(&tfm);

        let manual_tfm_extent =
            Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 1, 3].into());
        let mut manual_tfm_map = VecLatticeMap::fill(manual_tfm_extent, (0, 0, 0));
        for p in &manual_tfm_extent {
            // Have to do the inverse rotation here.
            *manual_tfm_map.get_world_ref_mut(&p) = (-p.y, p.x, p.z);
        }

        assert_eq!(tfm_map, manual_tfm_map);
    }

    #[test]
    fn test_rotationally_symmetric_map_serializes_equivalently() {
        let orig_extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 2, 2].into());
        let bottom_extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), [2, 2, 1].into());
        let top_extent = Extent::from_min_and_local_supremum([0, 0, 1].into(), [2, 2, 1].into());
        let mut orig_map = VecLatticeMap::<_, YLevelsIndexer>::fill(orig_extent, 0);
        for p in &bottom_extent {
            *orig_map.get_world_ref_mut(&p) = 1;
        }
        for p in &top_extent {
            *orig_map.get_world_ref_mut(&p) = 2;
        }

        let orig_serial = orig_map.linearize_extent(&orig_extent);

        let matrix = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]];
        let tfm = Transform { matrix };

        // Apply successive rotations and make sure the serialization is always the same.
        let mut prev_map = orig_map;
        for _ in 0..4 {
            let tfm_map = prev_map.apply_octahedral_transform(&tfm);
            let tfm_serial = tfm_map.linearize_extent(&tfm_map.get_extent());
            assert_eq!(tfm_serial, orig_serial);
            prev_map = tfm_map;
        }
    }

    #[test]
    fn test_into_iter() {
        let extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), [1, 1, 1].into());
        let map = VecLatticeMap::<_, YLevelsIndexer>::fill(extent, 0);

        for (p, v) in LatticeMapIter(&map) {
            assert_eq!(p, [0, 0, 0].into());
            assert_eq!(*v, 0);
        }
    }
}

#[cfg(test)]
mod compression_tests {
    use super::*;

    use std::io::Write;

    const LZ4_LEVEL: u32 = 1;
    const MAP_SIZE: i32 = 128;

    #[test]
    fn portable_compress_and_decompress_benchmark() {
        let extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), [MAP_SIZE; 3].into());
        let map = VecLatticeMap::<_, YLevelsIndexer>::fill(extent, 0);

        let start = std::time::Instant::now();
        let compressed_map = map.compress(PortableLZ4 { level: LZ4_LEVEL });
        let elapsed_micros = start.elapsed().as_micros();
        std::io::stdout()
            .write(format!("portable compressing map took {} micros\n", elapsed_micros).as_bytes())
            .unwrap();

        let start = std::time::Instant::now();
        let decompressed_map = compressed_map.decompress();
        let elapsed_micros = start.elapsed().as_micros();
        std::io::stdout()
            .write(
                format!(
                    "portable decompressing map took {} micros\n",
                    elapsed_micros
                )
                .as_bytes(),
            )
            .unwrap();

        assert_eq!(map, decompressed_map);
    }

    #[test]
    fn fast_compress_and_decompress_benchmark() {
        let extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), [MAP_SIZE; 3].into());
        let map = VecLatticeMap::<_, YLevelsIndexer>::fill(extent, 0);

        let start = std::time::Instant::now();
        let compressed_map = map.compress(FastLZ4 { level: LZ4_LEVEL });
        let elapsed_micros = start.elapsed().as_micros();
        std::io::stdout()
            .write(
                format!(
                    "unportable compressing map took {} micros\n",
                    elapsed_micros
                )
                .as_bytes(),
            )
            .unwrap();

        let start = std::time::Instant::now();
        let decompressed_map = compressed_map.decompress();
        let elapsed_micros = start.elapsed().as_micros();
        std::io::stdout()
            .write(
                format!(
                    "unportable decompressing map took {} micros\n",
                    elapsed_micros
                )
                .as_bytes(),
            )
            .unwrap();

        assert_eq!(map, decompressed_map);
    }
}
