from astropy.io import fits
from astropy.utils.console import ProgressBar
from astrodendro import Dendrogram, pp_catalog
from astropy import units as u
import radio_beam
from astropy import wcs
import numpy as np
from matplotlib import pyplot as plt
from func import savereg, grabfileinfo
import argparse
import warnings
warnings.filterwarnings('ignore')

    
def detect(infile, region, band, min_value=0.000325, min_delta=0.0005525, min_npix=7.5, plot=False, verbose=True):
    
    outfile = 'region{}_band{}_val{:.5g}_delt{:.5g}_pix{}'.format(region, band, min_value, min_delta, min_npix)
    contfile = fits.open(infile)                        # load in fits image
    da = contfile[0].data.squeeze()                     # get rid of extra axes

    mywcs = wcs.WCS(contfile[0].header).celestial       # set up world coordinate system, ditch extra dimensions
    beam = radio_beam.Beam.from_fits_header(contfile[0].header)

    d = Dendrogram.compute(da, min_value=min_value, min_delta=min_delta, min_npix=min_npix, wcs=mywcs, verbose=verbose)
    pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg

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
    leaves = []
    for i in range(len(d.leaves)):
        leaves.append(d.leaves[i].idx)
    cat = cat[leaves]                                   # Use only leaves
    cat['_idx'] = range(len(cat['_idx']))               # Reassign star ids to be continuous
    
    # Use FWHM for ellipse dimensions instead of sigma
    cat['major_sigma'] = cat['major_sigma']*np.sqrt(8*np.log(2))
    cat['minor_sigma'] = cat['minor_sigma']*np.sqrt(8*np.log(2))
    cat.rename_column('major_sigma', 'major_fwhm')
    cat.rename_column('minor_sigma', 'minor_fwhm')
    cat.rename_column('flux', 'dend_flux_band{}'.format(band))
    
    # Rename _idx to include the band number in the hundreds digit
    for i in range(len(cat)):
        cat['_idx'][i] = int('{}{:02d}'.format(band, cat['_idx'][i]))
    
    # Output the catalog and region files
    cat.write('./cat/cat_'+outfile+'.dat', format='ascii')
    savereg(cat, './reg/reg_'+outfile+'.reg')

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


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Detect sources in an image file using dendrograms')
    parser.add_argument('region', metavar='region', type=str, help='name of the region as listed in "imgfileinfo.dat"')
    parser.add_argument('band', metavar='band', type=int, help='integer representing the ALMA band of observation')
    args = parser.parse_args()
    region = str(args.region)
    band = args.band
    
    infilename, _, _, min_value, delta_fraction, min_npix = grabfileinfo(region, band)

    print("Min value: ", min_value)
    print("Min delta: {:.5g}".format(min_value*delta_fraction))
    print("Min npix: ", min_npix, '\n')

    detect(infilename, region, band, min_value=min_value, min_delta=min_value*delta_fraction, min_npix=min_npix, verbose=True)
    
