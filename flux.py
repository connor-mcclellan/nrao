from glob import glob
from func import grabfileinfo
from astropy.io import fits
from astropy.table import Table, MaskedColumn
import regions
import astropy.units as u
from astropy import coordinates
from astropy.nddata.utils import Cutout2D, NoOverlapError
from func import rms, mask, apsum
import argparse
import numpy as np
from astropy import wcs
import radio_beam
from astropy.utils.console import ProgressBar

def flux(region):

    # import the catalog file, get names of bands
    filename = glob('./cat/mastercat_region{}*'.format(region))[0]
    catalog = Table(Table.read(filename, format='ascii'), masked=True)
    
    bands = np.array(filename.split('bands_')[1].split('.dat')[0].split('_'), dtype=int)
    n_bands = len(bands)
    n_rows = len(catalog)
    
    ellipse_npix_col = MaskedColumn(length=len(catalog), name='ellipse_npix', mask=True)
    circ1_npix_col = MaskedColumn(length=len(catalog), name='circ1_npix', mask=True)
    circ2_npix_col = MaskedColumn(length=len(catalog), name='circ2_npix', mask=True)
    circ3_npix_col = MaskedColumn(length=len(catalog), name='circ3_npix', mask=True)
    
    
    for i in range(n_bands):
        
        # Load image file for this band
        print("\nLoading image file for region {} in band {} (Image {}/{})".format(region, bands[i], i+1, n_bands))
        imfile = grabfileinfo(region, bands[i])[0]
        contfile = fits.open(imfile)
        data = contfile[0].data.squeeze()
        
        # Set up wcs, beam, and pixel scale for this image
        mywcs = wcs.WCS(contfile[0].header).celestial
        beam = radio_beam.Beam.from_fits_header(contfile[0].header)
        pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
        ppbeam = (beam.sr/(pixel_scale**2)).decompose().value
        print('ppbeam: ', ppbeam)
        data = data/ppbeam
        
        # Set up columns for each aperture
        peak_flux_col = MaskedColumn(length=len(catalog), name='peak_flux_band{}'.format(bands[i]), mask=True)
        
        ellipse_flux_col = MaskedColumn(length=len(catalog), name='ellipse_flux_band{}'.format(bands[i]), mask=True)
        circ1_flux_col = MaskedColumn(length=len(catalog), name='circ1_flux_band{}'.format(bands[i]), mask=True)
        circ2_flux_col = MaskedColumn(length=len(catalog), name='circ2_flux_band{}'.format(bands[i]), mask=True)
        circ3_flux_col = MaskedColumn(length=len(catalog), name='circ3_flux_band{}'.format(bands[i]), mask=True)
        
        ellipse_rms_col = MaskedColumn(length=len(catalog), name='ellipse_rms_band{}'.format(bands[i]), mask=True)
        circ1_rms_col = MaskedColumn(length=len(catalog), name='circ1_rms_band{}'.format(bands[i]), mask=True)
        circ2_rms_col = MaskedColumn(length=len(catalog), name='circ2_rms_band{}'.format(bands[i]), mask=True)
        circ3_rms_col = MaskedColumn(length=len(catalog), name='circ3_rms_band{}'.format(bands[i]), mask=True)
        
        circ1_r, circ2_r, circ3_r = 30, 50, 70
        
        print('Photometering sources')
        pb = ProgressBar(n_rows)
        
        # Iterate over sources, extracting ellipse parameters
        for j in range(n_rows):
            source = catalog[j]
            x_cen = source['x_cen']*u.deg
            y_cen = source['y_cen']*u.deg
            major = source['major_fwhm']*u.deg
            minor = source['minor_fwhm']*u.deg
            pa = source['position_angle']*u.deg
            
            # Convert to pixel coordinates
            position = coordinates.SkyCoord(x_cen, y_cen, frame='icrs', unit=(u.deg,u.deg))
            pix_position = np.array(position.to_pixel(mywcs))
            pix_major = major/pixel_scale
            pix_minor = minor/pixel_scale
            
            # Create cutout
            size = np.max([circ3_r*pixel_scale.value, major.value])*2.2*u.deg
            try:
                cutout = Cutout2D(data, position, size, mywcs, mode='partial')
            except NoOverlapError:
                pb.update()
                continue
            cutout_center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])
            
            # create all aperture shapes
            ellipse_reg = regions.EllipsePixelRegion(cutout_center, pix_major*2., pix_minor*2., angle=pa)
            circ1_reg = regions.CirclePixelRegion(cutout_center, circ1_r)
            circ2_reg = regions.CirclePixelRegion(cutout_center, circ2_r)
            circ3_reg = regions.CirclePixelRegion(cutout_center, circ3_r)
            
            # get flux information from regions
            ellipse_flux, ellipse_rms, peak_flux, ellipse_mask, ellipse_npix = apsum(ellipse_reg, cutout)
            circ1_flux, circ1_rms, _, circ1_mask, circ1_npix = apsum(circ1_reg, cutout)
            circ2_flux, circ2_rms, _, circ2_mask, circ2_npix = apsum(circ2_reg, cutout)
            circ3_flux, circ3_rms, _, circ3_mask, circ3_npix = apsum(circ3_reg, cutout)
            
            # add fluxes to appropriate columns
            peak_flux_col[j] = peak_flux
            ellipse_flux_col[j], ellipse_rms_col[j] = ellipse_flux, ellipse_rms
            circ1_flux_col[j], circ1_rms_col[j] = circ1_flux, circ1_rms
            circ2_flux_col[j], circ2_rms_col[j] = circ2_flux, circ2_rms
            circ3_flux_col[j], circ3_rms_col[j] = circ3_flux, circ3_rms
            
            ellipse_npix_col[j] = ellipse_npix
            circ1_npix_col[j] = circ1_npix
            circ2_npix_col[j] = circ2_npix
            circ3_npix_col[j] = circ3_npix
            
            pb.update()
    
        # add columns to catalog
        catalog.add_columns([peak_flux_col,
                             ellipse_flux_col, ellipse_rms_col,
                             circ1_flux_col, circ1_rms_col,
                             circ2_flux_col, circ2_rms_col,
                             circ3_flux_col, circ3_rms_col,])
                             
    catalog.add_columns([ellipse_npix_col, circ1_npix_col, circ2_npix_col, circ3_npix_col])
    
    # save catalog
    catalog.write(filename.split('.dat')[0]+'_photometered.dat', format='ascii')
    print("\nMaster catalog saved as '{}'".format(filename.split('.dat')[0]+'_photometered.dat'))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform aperture photometry on all sources in a master catalog.')
    parser.add_argument('region', metavar='region', type=str, help='name of the region as listed in "imgfileinfo.dat"')
    args = parser.parse_args()
    region = str(args.region)

    flux(region)


