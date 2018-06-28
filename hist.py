from utils import grabfileinfo, filter_masked, grabbands
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import matplotlib.gridspec as gs
import argparse

def hist(region, bands, shapes):

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    filename = glob('./cat/mastercat_region{}*_photometered.dat'.format(region))[0]
    t = Table(Table.read(filename, format='ascii'), masked=True)
    
    t = filter_masked(t, shapes)
    
    flux_array = np.ndarray((len(bands), len(shapes))) # bands are across rows, shapes are across columns
    for i in range(flux_array.shape[0]):
        for j in range(flux_array.shape[1]):
            flux_array[i][j] = np.array(t['{}_flux_band{}'.format(shapes[i], bands[j])])
    
    n_images = len(shapes)
    xplots = int(np.around(np.sqrt(n_images)))
    yplots = xplots
    fig = plt.figure(figsize=(12, 12))
    for l, shape in enumerate(shapes):
        ax = fig.add_subplots(xplots, yplots, i+1)
        for k, band in enumerate(bands):
            ax.hist(flux_array[k][l], alpha=0.5, color=colors[l])
    
            
            
hist('w51e2', [3, 6], ['ellipse', 'circ1', 'circ2', 'circ3'])
