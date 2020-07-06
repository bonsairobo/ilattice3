use crate::{Extent, GetExtent, Indexer, Lattice};

use image::{Rgba, RgbaImage};

impl<I: Indexer> Into<RgbaImage> for &Lattice<Rgba<u8>, I> {
    fn into(self) -> RgbaImage {
        let extent = self.get_extent();
        let size = extent.get_local_supremum();
        assert_eq!(size.z, 1);
        assert!(size.x > 0);
        assert!(size.y > 0);
        let (width, height) = (size.x as u32, size.y as u32);

        let mut img = RgbaImage::new(width, height);
        for p in extent {
            *img.get_pixel_mut(p.x as u32, p.y as u32) = *self.get_local(&p);
        }

        img
    }
}

impl<I: Indexer> From<(&RgbaImage, I)> for Lattice<Rgba<u8>, I> {
    fn from((image, indexer): (&RgbaImage, I)) -> Self {
        let size = [image.width() as i32, image.height() as i32, 1].into();
        let extent = Extent::from_min_and_local_supremum([0, 0, 0].into(), size);
        let mut lattice = Lattice::fill_with_indexer(indexer, extent, Rgba([0; 4]));
        for (x, y, pixel) in image.enumerate_pixels() {
            let point = [x as i32, y as i32, 0].into();
            *lattice.get_mut_local(&point) = *pixel;
        }

        lattice
    }
}
