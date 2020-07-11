use crate::{prelude::*, Extent, Indexer, Point, VecLatticeMap};

use dot_vox::*;

impl<I: Indexer> Into<DotVoxData> for VecLatticeMap<VoxColor, I> {
    fn into(self: Self) -> DotVoxData {
        let size = *self.get_extent().get_local_supremum();
        // Voxel coordinates are limited to u8.
        assert!(size.x <= std::u8::MAX as i32);
        assert!(size.y <= std::u8::MAX as i32);
        assert!(size.z <= std::u8::MAX as i32);

        let size = dot_vox::Size {
            x: size.x as u32,
            y: size.y as u32,
            z: size.z as u32,
        };

        let mut voxels = Vec::new();
        for p in self.get_extent() {
            let i = self.get_local(&p);
            if i != EMPTY_VOX_COLOR {
                assert!(i <= std::u8::MAX as VoxColor);
                let i = i as u8;
                voxels.push(dot_vox::Voxel {
                    x: p.x as u8,
                    y: p.y as u8,
                    z: p.z as u8,
                    i,
                });
            }
        }

        let model = dot_vox::Model { size, voxels };

        DotVoxData {
            version: 150,
            models: vec![model],
            palette: Vec::new(),
            materials: Vec::new(),
        }
    }
}

pub type VoxColor = u16;
pub const EMPTY_VOX_COLOR: VoxColor = std::u8::MAX as u16 + 1;

impl<I: Indexer> VecLatticeMap<VoxColor, I> {
    pub fn from_vox_with_indexer(indexer: I, data: &DotVoxData, model_index: usize) -> Self {
        let DotVoxData { models, .. } = data;
        let Model {
            size: Size { x, y, z },
            voxels,
        } = &models[model_index];
        let size = Point::new(*x as i32, *y as i32, *z as i32);
        let extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), size);
        let mut map = VecLatticeMap::fill_with_indexer(indexer, extent, EMPTY_VOX_COLOR);
        for Voxel { x, y, z, i } in voxels.into_iter() {
            let point = [*x as i32, *y as i32, *z as i32].into();
            *map.get_local_ref_mut(&point) = *i as VoxColor;
        }

        map
    }
}
