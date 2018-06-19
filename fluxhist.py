from astropy.table import Table
from matplotlib import pyplot as plt

def fluxhist(catfile):
    catalog = Table.read(catfile, format='ascii')
    circ_sums = catalog['circ_flux']
    dend_sums = catalog['dend_flux']
    
    #fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(12, 8))
    fig, ax1 = plt.subplots(1, 1)
    ax1.hist(circ_sums, bins=40, label='Circular Aperture Sums')
    ax1.set_ylabel('n')
    ax1.set_xlabel('Flux (Jy)')
    ax1.hist(dend_sums, bins=40, label='Dendrogram Contour Sums', alpha = 0.8)
    ax1.set_title('Flux Histogram')
    ax1.set_xlabel('Flux (Jy)')
    plt.legend()
    plt.show()
    
fluxhist('./cat/cat_val0.000325_delt0.0005525_pix7.5_filtered.dat')
