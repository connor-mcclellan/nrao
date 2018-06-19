import numpy as np
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
import pickle

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
            plt.imshow(image, origin='lower', cmap='gray', picker=True)
        else:
            plt.imshow(image, origin='lower', picker=True)
        plt.imshow(masks[i], origin='lower', cmap='gray', alpha=0.2)
        plt.text(0, 0, '{}  SN {:.1f}'.format(names[i], snr_vals[i]), fontsize=7, color='w')
        plt.xticks([])
        plt.yticks([])

def rms(x):
    return (np.absolute(np.mean(x**2) - (np.mean(x))**2))**0.5
    
def saveobj(obj, name):
    with open('./tests/'+name+'.pickle', 'wb') as output:
        pickle.dump(obj, output, protocol=pickle.HIGHEST_PROTOCOL)

def loadobj(name):
    filename = name+'.pickle'
    with open(filename, 'rb') as f:
        return pickle.load(f)
