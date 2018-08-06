import dendrocat
from astropy.io import fits
from astropy.table import Table
import argparse
from glob import glob
import time
import astropy.units as u
from astropy.utils.console import ProgressBar
import matplotlib.pyplot as plt
import regions
import numpy as np

def timecheck(action):
    start = time.time()
    action
    stop = time.time()
    return (stop-start)*u.second

def mass_photometry(fname, outfile):
    filenames = glob('/lustre/aoc/students/bmcclell/w51/'+fname)
    
    objs = []
    print('Loading files')
    pb = ProgressBar(len(filenames))
    for f in filenames:
        rs = dendrocat.RadioSource(fits.open(f))
        objs.append(rs)
        pb.update()
        
    #n = np.shape(objs[0].data)[0]
    #center = regions.PixCoord(n/2, n/2)
    #reg = regions.CirclePixelRegion(center, 3200)
    #mask = reg.to_mask()
    #img = mask.to_image((n, n)).astype(bool)
    #objs[0].data = np.where(img==True, objs[0].data, np.nan)
    
    # Debugging
    #plt.figure()
    #plt.imshow(objs[0].data)
    #plt.show()
    #objs[0].threshold = 4.5
    #objs[0].min_value = 1.1e-4
    #objs[0].min_delta = 1.2*objs[0].min_value
    #objs[0].to_catalog()
    #objs[0].autoreject()
    #objs[0].reject([44001, 44032])
    #dendrocat.utils.save_regions(objs[0].catalog, '/users/bmcclell/nrao/reg/test_mass_photometry.reg')
    #print('Autorejection complete')
    
    t = Table.read('/users/bmcclell/nrao/cat/w51IRS2_photometered.dat', format='ascii')
    #for col in t.colnames:
    #    if 'GHz' in col:
    #        t.remove_column(col)
    
    mc = dendrocat.MasterCatalog(*objs, catalog=t)
    print('\nMaster Catalog made')
    start = time.time()
    mc.photometer(dendrocat.ellipse)
    stop = time.time()
    print('Ellipse apertures photometered. Time: {} s'.format(stop-start))
    start = time.time()
    mc.photometer(dendrocat.annulus)
    stop = time.time()
    print('Annulus apertures photometered. Time: {} s'.format(stop-start))
    start = time.time()
    mc.catalog.write(outfile, format='ascii', overwrite=True)
    stop = time.time()
    print('Catalog written. Time: {} s'.format(stop-start))

parser = argparse.ArgumentParser(description='Photometer a large number of files at once.')
parser.add_argument('filename', metavar='filename', type=str, help='file names to grab using glob')
parser.add_argument('outfile', metavar='outfile', type=str, help='full file path to save output catalog to')
args = parser.parse_args()

fname = args.filename
outfile = args.outfile

mass_photometry(fname, outfile)
