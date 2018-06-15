from astropy.io import fits
from astropy.utils.console import ProgressBar
import os
from astrodendro import Dendrogram, pp_catalog
from astropy import units as u
from astropy import coordinates
from astropy.nddata.utils import Cutout2D
import radio_beam
from astropy import wcs
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gs
from copy import deepcopy
import warnings
warnings.filterwarnings('ignore')


def rms(x):
    return (np.absolute(np.mean(x**2) - (np.mean(x))**2))**0.5

def plot_grid(datacube, masks, rejects, snr_vals, names):
    n_images = len(datacube)
    xplots = int(np.around(np.sqrt(n_images)))
    yplots = xplots + 1
    gs1 = gs.GridSpec(yplots, xplots, wspace=0.0, hspace=0.0, top=1.-0.5/(xplots+1), bottom=0.5/(xplots+1), left=0.5/(yplots+1), right=1-0.5/(yplots+1))
    plt.figure(figsize=(9.5, 10))
    for i in range(n_images):
        image = datacube[i]
        plt.subplot(gs1[i])
        if rejects[i]:
            plt.imshow(image, origin='lower', cmap='gray')
        else:
            plt.imshow(image, origin='lower')
        plt.imshow(masks[i], origin='lower', cmap='gray', alpha=0.2)
        plt.text(0, 0, '{}  SN {:.1f}'.format(names[i], snr_vals[i]), fontsize=7, color='w')
        plt.xticks([])
        plt.yticks([])
    
def contour(infile, min_value=0.00035, min_delta=0.000525, min_npix=10, plot=True, verbose=True):
    
    outfile = 'dend_val{:.5g}_delt{:.5g}_pix{}'.format(min_value, min_delta, min_npix)
    contfile = fits.open(infile)                        # load in fits image
    da = contfile[0].data.squeeze()                     # get rid of extra axes

    mywcs = wcs.WCS(contfile[0].header).celestial       # set up world coordinate system, ditch extra dimensions
    beam = radio_beam.Beam.from_fits_header(contfile[0].header)

    d = Dendrogram.compute(da, min_value=min_value, min_delta=min_delta, min_npix=min_npix, wcs=mywcs, verbose=verbose)
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

    cat = pp_catalog(d, metadata)                       # set up position-position catalog
    if verbose:
        cat.pprint(show_unit=True, max_lines=10)        # display to check values
    
    with open('./reg/reg_'+outfile+'.reg', 'w') as fh:  # write catalog information to region file
	    fh.write("fk5\n")
	    for row in cat:
	        fh.write("ellipse({x_cen}, {y_cen}, {major_sigma}, {minor_sigma}, {position_angle}) # text={{{_idx}}}\n".format(**dict(zip(row.colnames, row))))

    if plot:                                            # create PDF plots of contour regions, if enabled
        ax = plt.gca()
        ax.cla()
        plt.imshow(da, cmap='gray_r', interpolation='none', origin='lower',
                  vmax=0.01, vmin=-0.001)
        pltr = d.plotter()
        
        if verbose:
            print("Plotting contours to PDF...")
            pb = ProgressBar(len(d.leaves))
            
        for struct in d.leaves:                         # iterate over each of the leaf structures
            pltr.plot_contour(ax, structure=struct, colors=['r'],
                              linewidths=[0.9], zorder=5)
            if struct.parent:
                while struct.parent:
                    struct = struct.parent
                pltr.plot_contour(ax, structure=struct, colors=[(0,1,0,1)],
                                  linewidths=[0.5])
            if verbose:
                pb.update()
                
        cntr = plt.gca().collections

        plt.setp([x for x in cntr if x.get_color()[0,0] == 1], linewidth=0.25)
        plt.setp([x for x in cntr if x.get_color()[0,1] == 1], linewidth=0.25)
        plt.savefig('./contour/contour_'+outfile+'.pdf')
        plt.axis((1125.4006254228616, 1670.3650637799306,
                 1291.6829155596627, 1871.8063499397681))
        plt.setp([x for x in cntr if x.get_color()[0,0] == 1], linewidth=0.75) # Red
        plt.setp([x for x in cntr if x.get_color()[0,1] == 1], linewidth=0.5) # Green
        plt.savefig('./contour/contour_'+outfile+'zoom.pdf')
        

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
    
    print('Calculating RMS values within aperture annuli...')
    pb = ProgressBar(len(rows))
    
    data_cube = []
    masks = []
    rejects = []
    snr_vals = []
    rejection_threshold = 7.
    
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
    plt.suptitle('min_value={}, min_delta={}, min_npix={}, rejection_threshold={}'.format(min_value, min_delta, min_npix, rejection_threshold))
    plt.show()

