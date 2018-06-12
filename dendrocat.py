from astropy.io import fits
from astropy.utils.console import ProgressBar
import os
from astrodendro import Dendrogram, pp_catalog
from astropy import units as u
import radio_beam
from astropy import wcs
import numpy as np
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings('ignore')


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


