import numpy as np
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
import astropy.units as u
import pickle

def plot_grid(datacube, masks, rejects, snr_vals, names):   #change this to use SNR vals and names from catalog column?
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
        
def convolve(major1, minor1, pa1, major2, minor2, pa2):
    alpha = ((major1 * np.cos(pa1))**2 + (minor1 * np.sin(pa1))**2 +
             (major2 * np.cos(pa2))**2 + (minor2 * np.sin(pa2))**2)
    beta = ((major1 * np.sin(pa1))**2 + (minor1 * np.cos(pa1))**2 +
            (major2 * np.sin(pa2))**2 + (minor2 * np.cos(pa2))**2)
    gamma = (2 * ((minor1**2 - major1**2) * np.sin(pa1) * np.cos(pa1) +
                  (minor2**2 - major2**2) * np.sin(pa2) * np.cos(pa2)))
    s = alpha + beta
    t = np.sqrt((alpha - beta)**2 + gamma**2)
    new_major = np.sqrt(0.5 * (s + t))
    new_minor = np.sqrt(0.5 * (s - t))
    if np.isclose(((abs(gamma) + abs(alpha - beta))**0.5).to(u.arcsec).value, 1e-7):
        new_pa = 0.0 * u.deg
    else:
        new_pa = 0.5 * np.arctan2(-1. * gamma, alpha - beta)
    return new_major, new_minor, new_pa

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
