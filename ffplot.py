from utils import grabcatname
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from glob import glob

# work this in as an imgfileinfo.dat column
ppbeam1 = 125.0926865069147 # band 3
ppbeam2 = 101.7249966753436 # band 6

# work this in as an imgfileinfo.dat column
nu6 = 226092028953.4
nu3 = 92982346121.92

# Dendrogram (abandoned code, rework in if shape=='dendrogram')
# gets rows only where dendrogram entries are unmasked
#index = np.intersect1d(np.where(t['detected_band6']==1), np.where(t['detected_band3']==1))

namedict = {
    'ellipse':'Elliptical',
    'circ1':'Small Circular',
    'circ2':'Medium Circular',
    'circ3':'Large Circular'            
}

def specindex(nu1, nu2, f2, alpha):
    return f2*(nu1/nu2)**alpha

# Tomorrow: make this so you can plot each set of aperture sums (differing in shape) on the same plot in different colors
def ffplot(region, shape, band1, band2, nu1, nu2, log=True):
    filename = glob('./cat/mastercat_region{}*_photometered.dat'.format(region))[0]
    t = Table(Table.read(filename, format='ascii'), masked=True)
    
    # this needs a rework -- test only the columns the user needs data from
    cols = ['peak_flux_band3', 'ellipse_flux_band3',  'ellipse_rms_band3', 'circ1_flux_band3', 'circ1_rms_band3', 'circ2_flux_band3', 'circ2_rms_band3', 'circ3_flux_band3', 'circ3_rms_band3', 'peak_flux_band6', 'ellipse_flux_band6', 'ellipse_rms_band6', 'circ1_flux_band6', 'circ1_rms_band6', 'circ2_flux_band6', 'circ2_rms_band6', 'circ3_flux_band6', 'circ3_rms_band6']
    index = list(set(range(len(t)))^set(np.nonzero(t.mask[cols])[0]))
    t = t[index]
    
    flux_band1 = t[shape+'_flux_band'+str(band1)]
    flux_band2 = t[shape+'_flux_band'+str(band2)]
    npix = t[shape+'_npix']
    marker_labels = t['_idx']

    xfluxes = np.linspace(np.min(flux_band1), np.max(flux_band1), 100)
    yfluxes2 = specindex(nu2, nu1, xfluxes, 2)
    yfluxes3 = specindex(nu2, nu1, xfluxes, 3)

    plt.figure()
    plt.errorbar(flux_band1, flux_band2, xerr=flux_band1/np.sqrt(npix/ppbeam1), yerr=flux_band2/np.sqrt(npix/ppbeam2), fmt='o', ms=2, color='k', elinewidth=0.5, label=None)
    plt.plot(xfluxes, yfluxes2, label='Spectral Index = 2')
    plt.plot(xfluxes, yfluxes3, label='Spectral Index = 3')

    for i, label in enumerate(marker_labels):
        plt.annotate(label, (flux_band1[i], flux_band2[i]), size=8)

    if log:
        plt.xlabel('Log Band {} Flux'.format(band1))
        plt.ylabel('Log Band {} Flux'.format(band2))
        plt.xscale('log')
        plt.yscale('log')
    else:
        plt.xlabel('Band {} Flux'.format(band1))
        plt.ylabel('Band {} Flux'.format(band2))
        
    plt.title('Flux v. Flux in Bands {} and {} with {} Apertures'.format(band1, band2, namedict[shape]))
    plt.legend()
    
if __name__ == '__main__':
    region = 'w51e2'
    ffplot(region, 'circ3', 3, 6, nu3, nu6)
    plt.show()
