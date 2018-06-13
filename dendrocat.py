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
import warnings
warnings.filterwarnings('ignore')


def rms(x):
    return (np.absolute(np.mean(x**2) - (np.mean(x))**2))**0.5


def contour(infile, min_value=0.0001, min_delta=0.0002, min_npix=10, plot=True, verbose=True):
    
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
        

def bgrms(infile, regfile):
    
    contfile = fits.open(infile)                        
    data = contfile[0].data.squeeze()                   
    mywcs = wcs.WCS(contfile[0].header).celestial
    outfile = os.path.basename(regfile)[4:-4]
    rows = np.loadtxt(regfile, dtype=str, skiprows=1, delimiter=' # text={*}')
    
    if os.path.isfile('./reg/reg_'+outfile+'_annulus.reg'):
        os.remove('./reg/reg_'+outfile+'_annulus.reg')
    
    print('Calculating RMS values within aperture annuli...')
    pb = ProgressBar(len(rows))
    
    for i in range(len(rows)):
        s = rows[i][8:-2]                                   # trim to just numbers for all the ellipses
        coords = np.asarray(s.split(', '), dtype=float)     # split string into a list of ellipse coordinates
        
        # ELLIPSE PARAMETERS
        x_cen = coords[0] #* u.deg
        y_cen = coords[1] #* u.deg
        major_sigma = coords[2] #* u.deg
        minor_sigma = coords[3] #* u.deg
        position_angle = coords[4] #* u.deg
        
        # Convert ellipse parameters to pixel values        
        x_pix, y_pix = np.array(mywcs.wcs_world2pix(x_cen, y_cen, 1))
        pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 #* u.deg / u.pix
        pix_major_axis = major_sigma/pixel_scale
        
        # Cutout section of the image we care about, to speed up computation time
        size = major_sigma*2.2*u.deg
        position = coordinates.SkyCoord(x_cen, y_cen, frame='fk5', unit=(u.deg,u.deg))
        cutout = Cutout2D(data, position, size, mywcs, mode='partial')
        
        yy,xx = np.indices(cutout.data.shape)
        xcp, ycp = cutout.wcs.wcs_world2pix(position.ra.deg, position.dec.deg, 0)
        rgrid = ((yy-ycp)**2+(xx-xcp)**2)**0.5
        
        # Measure background RMS within a circular annulus
        annulus_width = 15 #* u.pix
        n = rgrid.shape[0]
        mask = np.zeros((n, n), dtype=bool)
        y, x = np.ogrid[-y_pix:n-y_pix, -x_pix:n-x_pix]
        inner_mask = x**2. + y**2. <= pix_major_axis**2
        outer_mask = x**2. + y**2. <= (pix_major_axis+annulus_width)**2
        
        mask[outer_mask] = 1
        mask[inner_mask] = 0
        
        # Checking the positions of the masks
        #plt.imshow(mask, origin='lower', cmap='gray')
        #plt.show()
        
        pix_locs = np.where(mask == 1)
        pix_values = data[pix_locs[1], pix_locs[0]]
        bg_rms = rms(pix_values)
        
        # Add circular annulus coordinates to a new region file
        with open('./reg/reg_'+outfile+'_annulus.reg', 'a') as fh:  # write catalog information to region file
            if i == 0:
	            fh.write("fk5\n")
            fh.write("annulus({}, {}, {}, {}) # text={{RMS: {:.1E}}}\n".format(x_cen, y_cen, major_sigma, major_sigma+annulus_width*pixel_scale, bg_rms))
        pb.update()


