import numpy as np
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
import astropy.units as u
from astropy.table import Table
import pickle
from radio_beam import Beams

def grabfileinfo(region, band):
    """
    Searches the imgfileinfo.dat file for the image with the given region and band. Returns a 6-element tuple with structure (filename, region, band, min_value, delt_frac, min_npix).
    """
    f = Table.read('./imgfileinfo.dat', format='ascii')
    info = f[np.intersect1d(np.where(f['region'] == region), np.where(f['band'] == band))]
    return tuple(info[0])

def grabcatname(region, band, flag=''):
    _, region, band, val, delt_frac, pix = grabfileinfo(region, band)
    if flag:
        catfile = './cat/cat_region{}_band{}_val{}_delt{:.5g}_pix{}_{}.dat'.format(region, band, val, delt_frac*val, pix, flag)
    else:
        catfile = './cat/cat_region{}_band{}_val{}_delt{:.5g}_pix{}.dat'.format(region, band, val, delt_frac*val, pix)
    return catfile
    
def plot_grid(datacube, masks, rejects, snr_vals, names):
    n_images = len(datacube)
    xplots = int(np.around(np.sqrt(n_images)))
    yplots = xplots + 1
    gs1 = gs.GridSpec(yplots, xplots, wspace=0.0, hspace=0.0, top=1.-0.5/(xplots+1), bottom=0.5/(xplots+1), left=0.5/(yplots+1), right=1-0.5/(yplots+1))
    plt.figure(figsize=(9.5, 10))
    for i in range(n_images):
        image = datacube[i]
        plt.subplot(gs1[i])
        if rejects[i]:
            plt.imshow(image, origin='lower', cmap='gray')
        else:
            plt.imshow(image, origin='lower')
        plt.imshow(masks[i], origin='lower', cmap='gray', alpha=0.2)
        plt.text(0, 0, '{}  SN {:.1f}'.format(names[i], snr_vals[i]), fontsize=7, color='w')
        plt.xticks([])
        plt.yticks([])
        
def commonbeam(major1, minor1, pa1, major2, minor2, pa2):
    """
    Create a smallest bounding ellipse around two other ellipses. Give ellipse dimensions as astropy units quantities.
    """
    somebeams = Beams([major1.to(u.arcsec), major2.to(u.arcsec)] * u.arcsec, [minor1.to(u.arcsec), minor2.to(u.arcsec)] * u.arcsec, [pa1.to(u.deg), pa2.to(u.deg)] * u.deg)
    common = somebeams.common_beam()
    new_major = common._major
    new_minor = common._minor
    new_pa = common._pa
    return new_major.to(u.deg), new_minor.to(u.deg), new_pa

def apsum(region, cutout):
    reg_mask = mask(region, cutout)
    pixels = cutout.data[reg_mask.astype('bool')]
    ap_rms = rms(pixels)
    ap_sum = np.sum(pixels)
    peak_flux = np.max(pixels)
    
    return ap_sum, ap_rms, peak_flux, reg_mask, len(pixels)

def savereg(cat, filename):
    with open(filename, 'w') as fh:
        fh.write("icrs\n")
        for row in cat:
            fh.write("ellipse({x_cen}, {y_cen}, {major_fwhm}, {minor_fwhm}, {position_angle}) # text={{{_idx}}}\n".format(**dict(zip(row.colnames, row))))

def mask(reg, cutout):
    n = cutout.shape[0]
    mask = reg.to_mask(mode='center')
    return np.array(mask.to_image((n, n)), dtype='int')

def rms(x):
    return (np.absolute(np.mean(x**2) - (np.mean(x))**2))**0.5
    
def saveobj(obj, name):
    with open('./tests/'+name+'.pickle', 'wb') as output:
        pickle.dump(obj, output, protocol=pickle.HIGHEST_PROTOCOL)

def loadobj(name):
    filename = name+'.pickle'
    with open(filename, 'rb') as f:
        return pickle.load(f)
