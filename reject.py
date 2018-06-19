from astropy.io import fits
from astropy.utils.console import ProgressBar
import os
from astropy import units as u
from astropy import coordinates
from astropy.table import Table, Column
from astropy.nddata.utils import Cutout2D
from astropy import wcs
import numpy as np
from matplotlib import pyplot as plt
from copy import deepcopy
import radio_beam
from func import rms, plot_grid
import warnings
warnings.filterwarnings('ignore')

def reject(imfile, catfile, threshold):
    
    contfile = fits.open(imfile)
    data = contfile[0].data.squeeze()                   
    mywcs = wcs.WCS(contfile[0].header).celestial
    outfile = os.path.basename(catfile).split('cat_')[1].split('.dat')[0]
    catalog = Table.read(catfile, format='ascii')
    
    beam = radio_beam.Beam.from_fits_header(contfile[0].header)
    pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
    ppbeam = (beam.sr/(pixel_scale**2)).decompose().value
    
    min_value = outfile.split('val')[1].split('_delt')[0]
    min_delta = outfile.split('delt')[1].split('_pix')[0]
    min_npix = outfile.split('pix')[1]
    
    if os.path.isfile('./reg/reg_'+outfile+'_annulus.reg'):
        os.remove('./reg/reg_'+outfile+'_annulus.reg')
    if os.path.isfile('./reg/reg_'+outfile+'_filtered.reg'):
        os.remove('./reg/reg_'+outfile+'_filtered.reg')
    
    # Load in manually accepted sources
    overridden = []
    if os.path.isfile('./override/override_'+outfile+'.txt'):
        overridden = np.loadtxt('./override/override_'+outfile+'.txt')
    print("\nManually accepted sources: ", overridden)
    
    print('Calculating RMS values within aperture annuli...')
    pb = ProgressBar(len(catalog))
    
    data_cube = []
    masks = []
    rejects = []
    snr_vals = []
    circ_sums = []
    
    for i in range(len(catalog)):
        x_cen = catalog['x_cen'][i] * u.deg
        y_cen = catalog['y_cen'][i] * u.deg
        major_fwhm = catalog['major_fwhm'][i] * u.deg
        minor_fwhm = catalog['minor_fwhm'][i] * u.deg
        position_angle = catalog['position_angle'][i]
        dend_flux = catalog['dend_flux'][i]
        
        annulus_width = 15 #* u.pix
        center_distance = 10 #* u.pix
        
        # Convert ellipse parameters to pixel values        
        x_pix, y_pix = np.array(mywcs.wcs_world2pix(x_cen, y_cen, 1))
        pix_major_axis = major_fwhm/pixel_scale
        
        # Cutout section of the image we care about, to speed up computation time
        size = ((center_distance+annulus_width)*pixel_scale+major_fwhm)*2.2
        position = coordinates.SkyCoord(x_cen, y_cen, frame='icrs', unit=(u.deg,u.deg))
        cutout = Cutout2D(data, position, size, mywcs, mode='partial')
        xx, yy = cutout.center_cutout   # Start using cutout coordinates
        
        # Measure background RMS within a circular annulus
        n = cutout.shape[0]
        mask = np.zeros(cutout.shape, dtype=bool)
        y, x = np.ogrid[-yy:n-yy, -xx:n-xx]
        inner_mask = x**2. + y**2. <= (center_distance+pix_major_axis)**2
        outer_mask = x**2. + y**2. <= (center_distance+pix_major_axis+annulus_width)**2
        circ_mask = x**2. + y**2. <= pix_major_axis**2
        
        mask[outer_mask] = 1
        mask[inner_mask] = 0
        
        # Plots of annulus regions
        data_cube.append(cutout.data)
        graph_mask = deepcopy(mask)
        graph_mask[circ_mask] = 1
        masks.append(graph_mask)
        
        # Calculate the RMS within the annulus region
        bg_rms = rms(cutout.data[np.where(mask == 1)])
        
        # Find peak flux within center circle
        mask[:,:] = 0
        circ_mask = x**2. + y**2. <= pix_major_axis**2
        mask[circ_mask] = 1
        peak_flux = np.max(cutout.data[np.where(mask == 1)])
        flux_rms_ratio = peak_flux / bg_rms
        snr_vals.append(flux_rms_ratio)
        
        # Sum the flux within the circular aperture
        circ_sum = np.sum(cutout.data[np.where(mask == 1)])/ppbeam
        circ_sums.append(circ_sum)
        
        # Reject bad sources below some SNR threshold
        rejected = False
        if flux_rms_ratio <= threshold:
            if i not in overridden:
                rejected = True
        rejects.append(rejected)
        
        # Add non-rejected source ellipses to a new region file
        fname = './reg/reg_'+outfile+'_filtered.reg'
        with open(fname, 'a') as fh:  # write catalog information to region file
            if os.stat(fname).st_size == 0:
                fh.write("icrs\n")
            if not rejected:
                fh.write("ellipse({}, {}, {}, {}, {}) # text={{{}}}\n".format(x_cen, y_cen, major_fwhm, minor_fwhm, position_angle, i))
        pb.update()
        
    # Plot the grid of sources
    plot_grid(data_cube, masks, rejects, snr_vals, range(len(catalog)))
    plt.suptitle('min_value={}, min_delta={}, min_npix={}, threshold={:.4f}'.format(min_value, min_delta, min_npix, threshold))
    print('Input a comma-separated list of sources to manually accept, then close the plot window. ')
    plt.show()

    overrides = input("\nPress enter to confirm the above. ").split(', ')
    print(overrides)
    
    # Save the manually accepted sources
    fname = './override/override_'+outfile+'.txt'
    with open(fname, 'a') as fh:
        for num in overrides:
            fh.write('\n'+str(num))
    print("Manual overrides written to './override/override_"+outfile+".txt'. New overrides will take effect the next time the rejection script is run.")
    
    # Save the filtered catalog with new columns for circular aperture flux sum and SNR
    catalog.remove_rows(rejects)
    circ_sums = [s for s in circ_sums if not rejects[circ_sums.index(s)]]
    snr_vals = [s for s in snr_vals if not rejects[snr_vals.index(s)]]
    catalog['_idx'] = range(len(catalog['_idx']))               # Reassign star ids to be continuous
    catalog.add_column(Column(circ_sums), index=catalog.colnames.index('dend_flux')+1, name='circ_flux')
    catalog.add_column(Column(snr_vals), index=catalog.colnames.index('circ_flux')+1, name='snr')
    catalog.write('./cat/cat_'+outfile+'_filtered.dat', format='ascii')
    
# Execute the script
imfile = '/lustre/aoc/students/bmcclell/w51/W51e2_cont_briggsSC_tclean.image.fits.gz'
catfile = './cat/cat_val0.000325_delt0.0005525_pix7.5.dat'

reject(imfile, catfile, 6.)
