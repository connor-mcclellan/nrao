from astropy.io import fits
from astrodendro import Dendrogram, pp_catalog
from astropy import units as u
import radio_beam
from astropy import wcs
import numpy as np


image = fits.open('../w51e2_sci.spw0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19.mfs.I.manual.image.tt0.pbcor.fits')
da = image[0].data.squeeze()

mywcs = wcs.WCS(image[0].header).celestial
beam = radio_beam.Beam.from_fits_header(image[0].header)

d = Dendrogram.compute(da, min_value=0.0005, min_delta=0.0010, min_npix=10, wcs=mywcs, verbose=True)
pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
pixel_scale_as = pixel_scale.to(u.arcsec).value

metadata = {
	    'data_unit':u.Jy / u.beam,
            'spatial_scale':pixel_scale,
	    'beam_major':beam.major,
            'beam_minor':beam.minor,
	    'wavelength':9.298234612192E+10*u.Hz,
	    'velocity_scale':u.km/u.s,
	    'wcs':mywcs,
            }

cat = pp_catalog(d, metadata)
cat.pprint(show_unit=True, max_lines=10)

with open("cores.reg", 'w') as fh:
    fh.write("fk5\n")
    for row in cat:
        fh.write("ellipse({x_cen}, {y_cen}, {major_sigma}, "
                 "{minor_sigma}, {position_angle}) # text={{{_idx}}}\n"
.format(**dict(zip(row.colnames, row))))
