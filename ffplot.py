from utils import grabcatname
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from glob import glob

region = 'w51e2'
filename = glob('./cat/mastercat_region{}*_photometered.dat'.format(region))[0]

ppbeam3 = 125.0926865069147
ppbeam6 = 101.7249966753436

nu6 = 226092028953.4
nu3 = 92982346121.92

t = Table(Table.read(filename, format='ascii'), masked=True)

# Get sources only where both were detected in the original dendrogram
#index = np.intersect1d(np.where(t['detected_band6']==1), np.where(t['detected_band3']==1))

# Get sources where certain columns are all unmasked
cols = ['peak_flux_band3', 'ellipse_flux_band3',  'ellipse_rms_band3', 'circ1_flux_band3', 'circ1_rms_band3', 'circ2_flux_band3', 'circ2_rms_band3', 'circ3_flux_band3', 'circ3_rms_band3', 'peak_flux_band6', 'ellipse_flux_band6', 'ellipse_rms_band6', 'circ1_flux_band6', 'circ1_rms_band6', 'circ2_flux_band6', 'circ2_rms_band6', 'circ3_flux_band6', 'circ3_rms_band6']
index = list(set(range(len(t)))^set(np.nonzero(t.mask[cols])[0]))  # remove all indices where there are nonzero masks (i.e. masked=True) in any of these columns 
t = t[index]

flux_band6 = t['ellipse_flux_band6']
flux_band3 = t['ellipse_flux_band3']
npix = t['ellipse_npix']
marker_labels = t['_idx']

def specindex(nu1, nu2, f2, alpha):
    return f2*(nu1/nu2)**alpha

xfluxes = np.linspace(np.min(flux_band3), np.max(flux_band3), 100)    # these wil be along the band 3 axis
yfluxes2 = specindex(nu6, nu3, xfluxes, 2)
yfluxes3 = specindex(nu6, nu3, xfluxes, 3)

plt.figure()
plt.errorbar(flux_band3, flux_band6, xerr=flux_band3/np.sqrt(npix/ppbeam3), yerr=flux_band6/np.sqrt(npix/ppbeam6), fmt='o', ms=2, color='k', elinewidth=0.5, label=None)
plt.plot(xfluxes, yfluxes2, label='Spectral Index = 2')
plt.plot(xfluxes, yfluxes3, label='Spectral Index = 3')

for i, label in enumerate(marker_labels):
    plt.annotate(label, (flux_band3[i], flux_band6[i]), size=8)

plt.xlabel('Log Band 3 Flux')
plt.ylabel('Log Band 6 Flux')
plt.xscale('log')
plt.yscale('log')
plt.title('Flux v. Flux in Bands 3 and 6 with Elliptical Apertures')
plt.legend()
plt.show()
