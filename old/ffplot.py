from utils import grabfileinfo, filter_masked
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import argparse

# Dendrogram (abandoned code, rework in if shape=='dendrogram')
# gets rows only where dendrogram entries are unmasked
#index = np.intersect1d(np.where(t['detected_band6']==1), np.where(t['detected_band3']==1))

namedict = {
    'ellipse':'Elliptical',
    'circ1':'Small Circular',
    'circ2':'Medium Circular',
    'circ3':'Large Circular'            
}

def specindex(nu1, nu2, f1, alpha):
    return f1*(nu2/nu1)**(alpha)


def ffplot(region, shapes, band1, band2, log=True, label=False, grid=True, peak=False):
    """Make a flux v. flux plot.
    
    Parameters
    ----------
    region : str
        The region identifier for the images
    shapes : list (dtype=str)
        A list of aperture shapes to compare
    band1 : int
        The ALMA band of observation whose flux will be on the x-axis
    band2 : int
        The ALMA band of observation whose flux will be on the y-axis
    
    Keyword Arguments
    -----------
    log : bool (default True)
        Display the plot on a log scale
    label : bool (default False)
        Draw ID labels for each source
    grid : bool (default True)
        Display a separate plot for each aperture
    peak : bool (default False)
        Use peak flux within ellipse instead of any aperture sums
    """
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    filename = glob('./cat/mastercat_region{}*_photometered.dat'.format(region))[0]
    t = Table(Table.read(filename, format='ascii'), masked=True)
    
    # get nu1, nu2, ppbeam1, ppbeam2 from imgfileinfo.dat
    nu1, ppbeam1 = grabfileinfo(region, band1)[-2:]
    nu2, ppbeam2 = grabfileinfo(region, band2)[-2:]
    
    # Filter to rows where all values the user needs for plotting are unmasked
    t = filter_masked(t, shapes)
    
    # Grab the flux data
    flux_band1 = []
    flux_band2 = []
    npix = []
    
    if peak:
        flux_band1.append(t['peak_flux_band'+str(band1)])
        flux_band2.append(t['peak_flux_band'+str(band2)])
    else:
        for shape in shapes:
            flux_band1.append(t[shape+'_flux_band'+str(band1)])
            flux_band2.append(t[shape+'_flux_band'+str(band2)])
            npix.append(t[shape+'_npix'])
    marker_labels = t['_idx']
    
    specindex_xflux = np.linspace(np.min(flux_band1), np.max(flux_band1), 10)
    #specindex0_yflux = specindex(nu1, nu2, specindex_xflux, 0.5)
    specindex1_yflux = specindex(nu1, nu2, specindex_xflux, 1)
    specindex2_yflux = specindex(nu1, nu2, specindex_xflux, 2)
    specindex3_yflux = specindex(nu1, nu2, specindex_xflux, 3)

    # ------ PLOTTING ------ needs to be cleaned up
    if grid:
        n_images = len(shapes)
        xplots = int(np.around(np.sqrt(n_images)))
        yplots = xplots
        fig, axes = plt.subplots(ncols=yplots, nrows=xplots, figsize=(12, 12))
        for i in range(len(shapes)):
            ax = np.ndarray.flatten(axes)[i]
            ax.errorbar(flux_band1[i], flux_band2[i], xerr=flux_band1[i]/np.sqrt(npix[i]/ppbeam1), yerr=flux_band2[i]/np.sqrt(npix[i]/ppbeam2), fmt='o', ms=2, alpha=0.75, elinewidth=0.5, color=colors[i], label='{} Aperture Sums'.format(namedict[shapes[i]]))
            #ax.plot(specindex_xflux, specindex0_yflux, '--', color='pink', label='Spectral Index = 0.5')
            ax.plot(specindex_xflux, specindex1_yflux, '--', color='m', label='Spectral Index = 1')
            ax.plot(specindex_xflux, specindex2_yflux, '--', color='purple', label='Spectral Index = 2')
            ax.plot(specindex_xflux, specindex3_yflux, '--', color='k', label='Spectral Index = 3')
            ax.set_xticks([])
            ax.set_yticks([])
                 
            ax.set_xlim([.6*np.min(flux_band1), 1.4*np.max(flux_band1)])
            ax.set_ylim([.1*np.min(flux_band2), 1.9*np.max(flux_band2)])
            
            if label:
                for j, label in enumerate(marker_labels):
                    ax.annotate(label, (flux_band1[i][j], flux_band2[i][j]), size=8)
            if log:
                ax.set_xlabel('Log Band {} Flux'.format(band1))
                ax.set_ylabel('Log Band {} Flux'.format(band2))
                ax.set_xscale('log')
                ax.set_yscale('log')
            else:
                ax.set_xlabel('Band {} Flux'.format(band1))
                ax.set_ylabel('Band {} Flux'.format(band2))
            ax.legend()
    else:
        plt.figure()
        if peak:
            plt.errorbar(flux_band1, flux_band2, fmt='o', ms=2, alpha=0.75, elinewidth=0.5, label='Peak Flux')
        else:
            for i in range(len(shapes)):
                plt.errorbar(flux_band1[i], flux_band2[i], xerr=flux_band1[i]/np.sqrt(npix[i]/ppbeam1), yerr=flux_band2[i]/np.sqrt(npix[i]/ppbeam2), fmt='o', ms=2, alpha=0.75, elinewidth=0.5, label='{} Aperture Sums'.format(namedict[shapes[i]]))
        plt.plot(specindex_xflux, specindex1_yflux, '--', color='m', label='Spectral Index = 1')
        plt.plot(specindex_xflux, specindex2_yflux, '--', color='purple', label='Spectral Index = 2')
        plt.plot(specindex_xflux, specindex3_yflux, '--', color='k', label='Spectral Index = 3')
        
        plt.xlim([.6*np.min(flux_band1), 1.4*np.max(flux_band1)])
        plt.ylim([.1*np.min(flux_band2), 1.9*np.max(flux_band2)])
        
        if label:
            for i in range(len(shapes)):
                for j, label in enumerate(marker_labels):
                    plt.annotate(label, (flux_band1[i][j], flux_band2[i][j]), size=8)
        if peak:
            if log:
                plt.xlabel('Log Band {} Peak Flux'.format(band1))
                plt.ylabel('Log Band {} Peak Flux'.format(band2))
                plt.xscale('log')
                plt.yscale('log')
            else:
                plt.xlabel('Band {} Peak Flux'.format(band1))
                plt.ylabel('Band {} Peak Flux'.format(band2))
        else:
            if log:
                plt.xlabel('Log Band {} Flux'.format(band1))
                plt.ylabel('Log Band {} Flux'.format(band2))
                plt.xscale('log')
                plt.yscale('log')
            else:
                plt.xlabel('Band {} Flux'.format(band1))
                plt.ylabel('Band {} Flux'.format(band2))
        plt.legend()
    
    plt.suptitle('Flux v. Flux For Region {} In Bands {} and {}'.format(region, band1, band2))
    # ------------------
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Produce a flux flux plot with spectral indices for any number of aperture shapes')
    parser.add_argument('region', metavar='region', type=str, help='name of the region as listed in "imgfileinfo.dat"')
    parser.add_argument('bands', metavar='bands', type=int, nargs='+', help='integers representing the ALMA bands of observation')
    parser.add_argument('-l', '--log', action='store_true')
    parser.add_argument('-g', '--grid', action='store_true')
    parser.add_argument('-n', '--number', action='store_true')
    parser.add_argument('-p', '--peak', action='store_true')
    args = parser.parse_args()
    region = str(args.region)
    bands = sorted(args.bands)
    log = bool(args.log)
    grid = bool(args.grid)
    number = bool(args.number)
    peak = bool(args.peak)

    region = 'w51e2'
    #ffplot(region, ['ellipse', 'circ1', 'circ2', 'circ3'], bands[0], bands[1], label=number, grid=grid, log=log, peak=peak)
    ffplot(region, ['ellipse'], bands[0], bands[1], label=number, grid=grid, log=log, peak=peak)
    plt.show()
    
