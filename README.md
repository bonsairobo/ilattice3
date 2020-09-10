# ilattice3: Voxel Structures

`ilattice3` provides many useful types and functions for dealing with voxel
data. The core types are `Point`, `Extent`, and the various "lattice maps,"
which allow storing (or generating), data at each point in a 3D integer lattice.

Consider this example which samples a signed distance function at each point
and stores it in a compact vector structure.

```rust
fn sphere_sdf(p: &Point) -> f32 {
    let d = p.dot(p);

    (d - 10) as f32
}

let sample_extent = Extent::from_center_and_radius([0, 0, 0].into(), 20);
let sampled_sphere =
    VecLatticeMap::copy_from_map(&FnLatticeMap::new(sphere_sdf), sample_extent);
```

In addition to the `FnLatticeMap` and `VecLatticeMap`, there are:
* `ChunkedLatticeMap`: a sparse lattice map implemented as a hash map of `VecLatticeMap`s
* `ChunkedPaletteLatticeMap`: a `ChunkedLatticeMap` with an additional layer of palette compression

Other miscellaneous features:
* fast LZ4 compression of `VecLatticeMap`
* conversion to/from VOX format
* conversion to/from RgbaImage format

There are also companion modules:
* [ilattice3-mesh](https://github.com/bonsairobo/ilattice3-mesh): various meshing algorithms for voxel data
* [ilattice3-wfc](https://github.com/bonsairobo/ilattice3-wfc): the Wave Function Collapse algorithm for 2D and 3D lattices

TODO:
* `ChunkedLatticeMap` with variable level of detail per chunk
