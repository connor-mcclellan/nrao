from glob import glob
from astropy.io import fits
from astropy.table import Table, MaskedColumn
import regions
import astropy.units as u
from astropy import coordinates
from astropy.nddata.utils import Cutout2D, NoOverlapError
from utils import rms, mask, apsum, grabfileinfo, plot_grid
import argparse
import numpy as np
from astropy import wcs
import radio_beam
from astropy.utils.console import ProgressBar
import matplotlib.pyplot as plt

def flux(region):

    # import the catalog file, get names of bands
    filename = glob('./cat/mastercat_region{}*'.format(region))[0]
    catalog = Table(Table.read(filename, format='ascii'), masked=True)
    catalog.sort('_idx')
    
    bands = np.array(filename.split('bands_')[1].split('.dat')[0].split('_'), dtype=int)
    n_bands = len(bands)
    n_rows = len(catalog)
    
    ellipse_npix_col = MaskedColumn(length=len(catalog), name='ellipse_npix', mask=True)
    circ1_npix_col = MaskedColumn(length=len(catalog), name='circ1_npix', mask=True)
    circ2_npix_col = MaskedColumn(length=len(catalog), name='circ2_npix', mask=True)
    circ3_npix_col = MaskedColumn(length=len(catalog), name='circ3_npix', mask=True)
    
    n_rejected = 0
    
    for i in range(n_bands):
        
        band = bands[i]
        
        # Load image file for this band
        print("\nLoading image file for region {} in band {} (Image {}/{})".format(region, band, i+1, n_bands))
        imfile = grabfileinfo(region, band)[0]
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
        peak_flux_col = MaskedColumn(length=len(catalog), name='peak_flux_band{}'.format(band), mask=True)
        annulus_median_col = MaskedColumn(length=len(catalog), name='annulus_median_band{}'.format(band), mask=True)
        annulus_rms_col = MaskedColumn(length=len(catalog), name='annulus_rms_band{}'.format(band), mask=True)
        ellipse_flux_col = MaskedColumn(length=len(catalog), name='ellipse_flux_band{}'.format(band), mask=True)
        circ1_flux_col = MaskedColumn(length=len(catalog), name='circ1_flux_band{}'.format(band), mask=True)
        circ2_flux_col = MaskedColumn(length=len(catalog), name='circ2_flux_band{}'.format(band), mask=True)
        circ3_flux_col = MaskedColumn(length=len(catalog), name='circ3_flux_band{}'.format(band), mask=True)
        
        ellipse_rms_col = MaskedColumn(length=len(catalog), name='ellipse_rms_band{}'.format(band), mask=True)
        circ1_rms_col = MaskedColumn(length=len(catalog), name='circ1_rms_band{}'.format(band), mask=True)
        circ2_rms_col = MaskedColumn(length=len(catalog), name='circ2_rms_band{}'.format(band), mask=True)
        circ3_rms_col = MaskedColumn(length=len(catalog), name='circ3_rms_band{}'.format(band), mask=True)
        
        circ1_r, circ2_r, circ3_r = 5e-6*u.deg, 1e-5*u.deg, 1.5e-5*u.deg
        
        print('Photometering sources')
        pb = ProgressBar(len(catalog[np.where(catalog['rejected']==0)]))
        
        masks = []
        datacube = []
        rejects = []
        snr_vals = []
        names = []
        
        # Iterate over sources, extracting ellipse parameters
        for j in range(n_rows):
        
            if catalog['rejected'][j] == 1:
                continue
        
            source = catalog[j]
            x_cen = source['x_cen']*u.deg
            y_cen = source['y_cen']*u.deg
            major = source['major_fwhm']*u.deg
            minor = source['minor_fwhm']*u.deg
            pa = source['position_angle']*u.deg
            
            annulus_width = 1e-5*u.deg
            center_distance = 1e-5*u.deg
            
            # Convert to pixel coordinates
            position = coordinates.SkyCoord(x_cen, y_cen, frame='icrs', unit=(u.deg,u.deg))
            pix_position = np.array(position.to_pixel(mywcs))
            pix_major = major/pixel_scale
            pix_minor = minor/pixel_scale
            
            # Create cutout
            size = np.max([circ3_r.value, major.value+center_distance.value+annulus_width.value])*2.2*u.deg
            try:
                cutout = Cutout2D(data, position, size, mywcs, mode='partial')
            except NoOverlapError:
                catalog['rejected'][j] = 1
                pb.update()
                continue
            cutout_center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])
            datacube.append(cutout.data)
            
            # create all aperture shapes
            ellipse_reg = regions.EllipsePixelRegion(cutout_center, pix_major*2., pix_minor*2., angle=pa)
            circ1_reg = regions.CirclePixelRegion(cutout_center, circ1_r/pixel_scale)
            circ2_reg = regions.CirclePixelRegion(cutout_center, circ2_r/pixel_scale)
            circ3_reg = regions.CirclePixelRegion(cutout_center, circ3_r/pixel_scale)
            
            innerann_reg = regions.CirclePixelRegion(cutout_center, center_distance/pixel_scale+pix_major)
            outerann_reg = regions.CirclePixelRegion(cutout_center, center_distance/pixel_scale+pix_major+annulus_width/pixel_scale)
            
            annulus_mask = mask(outerann_reg, cutout) - mask(innerann_reg, cutout)
            
            # get flux information from regions
            ellipse_flux, ellipse_rms, peak_flux, ellipse_mask, ellipse_npix = apsum(ellipse_reg, cutout)
            circ1_flux, circ1_rms, _, circ1_mask, circ1_npix = apsum(circ1_reg, cutout)
            circ2_flux, circ2_rms, _, circ2_mask, circ2_npix = apsum(circ2_reg, cutout)
            circ3_flux, circ3_rms, _, circ3_mask, circ3_npix = apsum(circ3_reg, cutout)
            
            annulus_rms = rms(cutout.data[annulus_mask.astype('bool')])
            annulus_median = np.median(cutout.data[annulus_mask.astype('bool')])
            
            # Add grid plot mask to list
            masklist = [ellipse_mask, annulus_mask, circ1_mask, circ2_mask, circ3_mask]
            masks.append(masklist)
            
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
            
            annulus_median_col[j] = annulus_median
            annulus_rms_col[j] = annulus_rms
            
            catalog['snr_band'+str(band)][j] = peak_flux/annulus_rms
            snr_vals.append(peak_flux/annulus_rms)
            names.append(catalog['_idx'][j])
            
            # Secondary rejection
            rejected = 0
            lowest_flux = np.min([ellipse_flux, circ1_flux, circ2_flux, circ3_flux])
            #if lowest_flux <= annulus_median*ellipse_npix or lowest_flux < 0:
            if lowest_flux < 0:
                catalog['rejected'][j] = 1
                n_rejected += 1
                rejected = 1
            rejects.append(rejected)
            pb.update()
    
        # Plot the grid of sources
        plot_grid(datacube, masks, rejects, snr_vals, names)
        plt.suptitle('region={}, band={}'.format(region, band))  
        plt.show(block=False)
    
        # add columns to catalog
        catalog.add_columns([peak_flux_col,
                             ellipse_flux_col, ellipse_rms_col,
                             circ1_flux_col, circ1_rms_col,
                             circ2_flux_col, circ2_rms_col,
                             circ3_flux_col, circ3_rms_col,])
                             
    catalog.add_columns([ellipse_npix_col, circ1_npix_col, circ2_npix_col, circ3_npix_col, annulus_median_col, annulus_rms_col])
    print("\n{} sources flagged for secondary rejection".format(n_rejected))
    
    # save catalog
    catalog = catalog[sorted(catalog.colnames)]
    catalog.write(filename.split('.dat')[0]+'_photometered.dat', format='ascii')
    print("\nMaster catalog saved as '{}'".format(filename.split('.dat')[0]+'_photometered.dat'))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform aperture photometry on all sources in a master catalog.')
    parser.add_argument('region', metavar='region', type=str, help='name of the region as listed in "imgfileinfo.dat"')
    args = parser.parse_args()
    region = str(args.region)

    flux(region)
    input("Press enter to terminate.")


