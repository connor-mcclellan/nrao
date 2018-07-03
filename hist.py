from utils import grabfileinfo, filter_masked, grabbands
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import argparse

namedict = {
    'ellipse':'Elliptical',
    'circ1':'Small Circular',
    'circ2':'Medium Circular',
    'circ3':'Large Circular'            
}

def hist(region, bands, shapes):

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    filename = glob('./cat/mastercat_region{}*_photometered.dat'.format(region))[0]
    t = Table(Table.read(filename, format='ascii'), masked=True)
    
    t = filter_masked(t, shapes)
    
    flux_array = np.ndarray((len(bands), len(shapes), len(t))) # bands are across rows, shapes are across columns, sources are across depth
    for i in range(len(bands)):
        for j in range(len(shapes)):
            flux_array[i][j] = np.array(t['{}_flux_band{}'.format(shapes[j], bands[i])])
    
    n_images = len(shapes)
    xplots = int(np.around(np.sqrt(n_images)))
    yplots = xplots
    fig = plt.figure(figsize=(6, 6))
    
    for l, shape in enumerate(shapes):
        ax = fig.add_subplot(xplots, yplots, l+1)
        for k, band in enumerate(bands):
            n, bins, patches = ax.hist(flux_array[k][l], bins=np.logspace(-5, 0, 20), alpha=0.5, color=colors[k+3], label='Band {}'.format(band))
        ax.set_title('{} Apertures'.format(namedict[shapes[l]]))
        ax.set_ylim(bottom=0, top=12)
        ax.set_xlabel('log10(Flux) (Jy)')
        ax.set_xscale('log')
        ax.set_ylabel('n')
        ax.legend()
    plt.suptitle('Flux Histograms For Region {} In Bands {}'.format(region, bands))
    plt.subplots_adjust(top=0.897,
bottom=0.108,
left=0.052,
right=0.985,
hspace=0.409,
wspace=0.219)

if __name__ == '__main__':
    hist('w51e2', [3, 6], ['ellipse', 'circ1', 'circ2', 'circ3'])
    plt.show()
