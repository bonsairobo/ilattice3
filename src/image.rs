use crate::{prelude::*, Extent, Indexer, VecLatticeMap};

use image::{GenericImageView, ImageBuffer, Pixel};

impl<P, I> Into<ImageBuffer<P, Vec<<P as Pixel>::Subpixel>>> for &VecLatticeMap<P, I>
where
    P: Pixel + 'static,
    I: Indexer,
{
    fn into(self) -> ImageBuffer<P, Vec<<P as Pixel>::Subpixel>> {
        let extent = self.get_extent();
        let size = extent.get_local_supremum();
        assert_eq!(size.z, 1);
        assert!(size.x > 0);
        assert!(size.y > 0);
        let (width, height) = (size.x as u32, size.y as u32);

        let mut img = ImageBuffer::new(width, height);
        for p in extent {
            *img.get_pixel_mut(p.x as u32, p.y as u32) = self.get_local(&p);
        }

        img
    }
}

impl<Im, I> From<(&Im, I)> for VecLatticeMap<<Im as GenericImageView>::Pixel, I>
where
    Im: GenericImageView,
    I: Indexer,
{
    fn from((image, _indexer): (&Im, I)) -> Self {
        let size = [image.width() as i32, image.height() as i32, 1].into();
        let extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), size);
        let mut map = unsafe { VecLatticeMap::uninitialized(extent) };
        for (x, y, pixel) in image.pixels() {
            let point = [x as i32, y as i32, 0].into();
            *map.get_local_ref_mut(&point) = pixel;
        }

        map
    }
}
