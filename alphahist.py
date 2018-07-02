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

def specindex(nu1, nu2, f1, f2):
    return np.log(nu1/nu2)/np.log(f1/f2)

def hist(region, bands, shapes):

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    filename = glob('./cat/mastercat_region{}*_photometered.dat'.format(region))[0]
    t = Table(Table.read(filename, format='ascii'), masked=True)
    
    t = filter_masked(t, shapes)
    
    # NOT GENERALIZED
    nu3, ppbeam3 = grabfileinfo(region, bands[0])[-2:]
    nu6, ppbeam6 = grabfileinfo(region, bands[1])[-2:]
    
    flux_array = np.ndarray((len(bands), len(shapes), len(t))) # bands are across rows, shapes are across columns, sources are across depth
    for i in range(len(bands)):
        for j in range(len(shapes)):
            flux_array[i][j] = np.array(t['{}_flux_band{}'.format(shapes[j], bands[i])])
    
    fig = plt.figure(figsize=(12, 6))
    
    for l, shape in enumerate(shapes):
        ax = fig.add_subplot(1, 2, l+1)
        print(nu3, nu6, flux_array[0][l][0]/flux_array[1][l][0], specindex(nu3, nu6, flux_array[0][l][0], flux_array[1][l][0]))
        n, bins, patches = ax.hist(specindex(nu3, nu6, flux_array[0][l], flux_array[1][l]), bins=np.linspace(0, 4, 10), color=colors[0])
        ax.set_title('{} Apertures'.format(namedict[shapes[l]]))
        ax.set_xlabel('Alpha')
        #ax.set_xscale('log')
        ax.set_ylabel('n')
    plt.suptitle('Alpha Parameter For Sources In Region {} In Bands {}'.format(region, bands))
    plt.subplots_adjust(top=0.887,
bottom=0.088,
left=0.052,
right=0.985,
hspace=0.234,
wspace=0.129)
            
hist('w51e2', [3, 6], ['circ1', 'circ2'])
plt.show()
