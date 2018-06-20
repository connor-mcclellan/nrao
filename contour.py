from astropy.io import fits
from astropy.utils.console import ProgressBar
from astrodendro import Dendrogram, pp_catalog
from astropy import units as u
import radio_beam
from astropy import wcs
import numpy as np
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings('ignore')

    
def contour(infile, region, band, min_value=0.000325, min_delta=0.0005525, min_npix=7.5, plot=False, verbose=False):
    
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
    cat.rename_column('flux', 'dend_flux')
    
    # Output the catalog file
    cat.write('./cat/cat_'+outfile+'.dat', format='ascii')
    
    with open('./reg/reg_'+outfile+'.reg', 'w') as fh:  # write catalog information to region file
	    fh.write("icrs\n")
	    for row in cat:
	        fh.write("ellipse({x_cen}, {y_cen}, {major_fwhm}, {minor_fwhm}, {position_angle}) # text={{{_idx}}}\n".format(**dict(zip(row.colnames, row))))

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
        plt.savefig('./plot_contour/contour_'+outfile+'.pdf')
        plt.axis((1125.4006254228616, 1670.3650637799306,
                 1291.6829155596627, 1871.8063499397681))
        plt.setp([x for x in cntr if x.get_color()[0,0] == 1], linewidth=0.75) # Red
        plt.setp([x for x in cntr if x.get_color()[0,1] == 1], linewidth=0.5) # Green
        plt.savefig('./plot_contour/contour_'+outfile+'zoom.pdf')
        
#--- SET IMAGE INFORMATION ---
        
#infilename = '/lustre/aoc/students/bmcclell/w51/w51e2_sci.spw0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19.mfs.I.manual.image.tt0.pbcor.fits.gz'    # band 3
#band = 3
#region = w51e2

infilename = '/lustre/aoc/students/bmcclell/w51/W51e2_cont_briggsSC_tclean.image.fits.gz'   # band 6
band = 6
region = w51e2

#-----------------------------

#min_values = np.linspace(0.00015, 0.000245, 6)
min_values = np.array([0.000325])
min_deltas = min_values*1.7
min_npixs = [7.5]

print("Min values: ", min_values)
print("Min deltas: ", min_deltas)
print("Min npix: ", min_npixs, '\n')

print('Total task progress:')
pb = ProgressBar(len(min_values)*len(min_npixs))
print()

for i in range(len(min_values)):
    for j in range(len(min_npixs)):
        print("\nMin value: {:.8g}    Min delta: {:.8g}   Min npix: {}".format(min_values[i], min_deltas[i], min_npixs[j]))
        contour(infilename, region, band, min_value=min_values[i], min_delta=min_deltas[i], min_npix=min_npixs[j], verbose=False)
        
        # print('Total task progress:')
        pb.update()
    
