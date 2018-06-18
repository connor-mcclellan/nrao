from astropy.io import fits
from astropy.utils.console import ProgressBar
import os
from astropy import units as u
from astropy import coordinates
from astropy.nddata.utils import Cutout2D
from astropy import wcs
import numpy as np
from matplotlib import pyplot as plt
from copy import deepcopy
from func import rms, plot_grid
import warnings
warnings.filterwarnings('ignore')

def reject(infile, regfile):
    
    contfile = fits.open(infile)                        
    data = contfile[0].data.squeeze()                   
    mywcs = wcs.WCS(contfile[0].header).celestial
    outfile = os.path.basename(regfile).split('reg_')[1].split('.reg')[0]
    rows = np.loadtxt(regfile, dtype=str, skiprows=1, delimiter=' # text={*}')
    
    min_value = outfile.split('val')[1].split('_delt')[0]
    min_delta = outfile.split('delt')[1].split('_pix')[0]
    min_npix = outfile.split('pix')[1]
    
    if os.path.isfile('./reg/reg_'+outfile+'_annulus.reg'):
        os.remove('./reg/reg_'+outfile+'_annulus.reg')
    if os.path.isfile('./reg/reg_'+outfile+'_filtered.reg'):
        os.remove('./reg/reg_'+outfile+'_filtered.reg')
    
    # Load in manually accepted sources
    accepted = []
    if os.path.isfile('./override/'+outfile+'_override.txt'):
        accepted = np.loadtxt('./override/'+outfile+'_override.txt')
    print("\nManually accepted sources: ", accepted)
    
    print('Calculating RMS values within aperture annuli...')
    pb = ProgressBar(len(rows))
    
    data_cube = []
    masks = []
    rejects = []
    snr_vals = []
    rejection_threshold = 3.51e4*rms(data)      # This may need some adjustment. Lower RMS means less noise, and so more sources should be accepted since they're less likely to be noise.
    
    for i in range(len(rows)):
        s = rows[i].split('ellipse(')[1].split(') ')[0]                                   # trim to just numbers for all the ellipses
        coords = np.asarray(s.split(', '), dtype=float)     # split string into a list of ellipse coordinates
        
        # ELLIPSE PARAMETERS
        x_cen = coords[0] #* u.deg
        y_cen = coords[1] #* u.deg
        major_sigma = coords[2] #* u.deg
        minor_sigma = coords[3] #* u.deg
        position_angle = coords[4] #* u.deg
        
        annulus_width = 15 #* u.pix
        center_distance = 10 #* u.pix
        
        # Convert ellipse parameters to pixel values        
        x_pix, y_pix = np.array(mywcs.wcs_world2pix(x_cen, y_cen, 1))
        pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 #* u.deg / u.pix
        pix_major_axis = major_sigma/pixel_scale
        
        # Cutout section of the image we care about, to speed up computation time
        size = ((center_distance+annulus_width)*pixel_scale+major_sigma)*2.2*u.deg
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
        
        # Reject bad sources below some SNR threshold
        rejected = False
        if flux_rms_ratio <= rejection_threshold:     # How should we determine this threshold? How does it relate to a source being a <some number> sigma detection?
            if i not in accepted:
                rejected = True
        rejects.append(rejected)
        
        # Add circular annulus coordinates to a new region file
        with open('./reg/reg_'+outfile+'_annulus.reg', 'a') as fh:  # write catalog information to region file
            if i == 0:
	            fh.write("icrs\n")
            fh.write("annulus({}, {}, {}, {}) # text={{#{} SNR: {:.1f}}}\n".format(x_cen, y_cen, center_distance*pixel_scale+major_sigma, pixel_scale*(center_distance+annulus_width)+major_sigma, i, flux_rms_ratio))
            
        # Add non-rejected source ellipses to a new region file
        if not rejected:
            fname = './reg/reg_'+outfile+'_filtered.reg'
            with open(fname, 'a') as fh:  # write catalog information to region file
                if os.stat(fname).st_size == 0:
	                fh.write("icrs\n")
                fh.write("ellipse({}, {}, {}, {}, {}) # text={{{}}}\n".format(x_cen, y_cen, major_sigma, minor_sigma, position_angle, i))
            
        pb.update()
        
    # Plot the grid of sources
    plot_grid(data_cube, masks, rejects, snr_vals, range(len(rows)))
    plt.suptitle('min_value={}, min_delta={}, min_npix={}, rejection_threshold={:.4f}'.format(min_value, min_delta, min_npix, rejection_threshold))
    print('Input a comma-separated list of sources to manually accept, then close the plot window. ')
    plt.show()

    overrides = input("\nPress enter to confirm the above. ").split(', ')
    print(overrides)
    
    fname = './override/'+outfile+'_override.txt'
    with open(fname, 'a') as fh:
        for num in overrides:
            fh.write('\n'+str(num))
    print("Manual overrides written to './override/"+outfile+"_override.txt'. New overrides will take effect the next time the rejection script is run.")
    
    
# Execute the script
infile = '/lustre/aoc/students/bmcclell/w51/W51e2_cont_briggsSC_tclean.image.fits.gz'
regfile = './reg/reg_dend_val0.000325_delt0.0004875_pix7.5.reg'

reject(infile, regfile)
