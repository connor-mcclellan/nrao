from astropy.table import Table
from matplotlib import pyplot as plt

def fluxhist(catfile):
    catalog = Table.read(catfile, format='ascii')
    circ_sums = catalog['circ_flux']
    dend_sums = catalog['dend_flux']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    ax1.hist(circ_sums, bins=40)
    ax1.set_title('Circular Aperture Sums')
    ax1.set_ylabel('n')
    ax1.set_xlabel('Flux (Jy)')
    ax2.hist(dend_sums, bins=40)
    ax2.set_title('Dendrogram Contour Sums')
    ax2.set_xlabel('Flux (Jy)')
    plt.show()
    
fluxhist('./cat/cat_val0.000325_delt0.0005525_pix7.5_filtered.dat')
