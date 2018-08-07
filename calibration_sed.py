from astropy.table import Table
from astropy.utils.console import ProgressBar
import dendrocat
import argparse
import matplotlib.pyplot as plt

def calib_sed(catalogname, outpath):

    t = Table.read(catalogname, format='ascii')
    t_accepted = t[t['rejected']==0]
    
    for i in range(len(t_accepted)):
        print(i)
        dendrocat.utils.plot_sed(t_accepted[i], t_accepted, alphas=[1,2,3])
        plt.show()
        
parser = argparse.ArgumentParser(description='Plot an SED for calibration data.')
parser.add_argument('catname', metavar='catname', type=str, help='filename to the photometry catalog')
parser.add_argument('outpath', metavar='outpath', type=str, help='file path to save output plots to')
args = parser.parse_args()

fname = args.catname
outfile = args.outpath

calib_sed(fname, outfile)
        
