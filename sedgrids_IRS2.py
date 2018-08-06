from astropy.table import Table
from astropy.io import fits
import dendrocat
import astropy.units as u 
import numpy as np
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
from collections import OrderedDict


irs2 = Table.read('/users/bmcclell/nrao/cat/w51IRS2_photometered.dat', format='ascii')

n1 = dendrocat.RadioSource(fits.open('/users/bmcclell/Data/irs2/w51n_band6_226GHz_adjustedRADEC.fits'))
n2 = dendrocat.RadioSource(fits.open('/users/bmcclell/Data/irs2/w51north_sci.spw0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19.mfs.I.manual.image.tt0.pbcor.fits'))
n3 = dendrocat.RadioSource(fits.open('/users/bmcclell/Data/irs2/w51n_bandQ_45GHz_adjustedRADEC.fits'))

n3.nu = 45.0*u.GHz
n3.freq_id = '45.0GHz'
n3.set_metadata()

mc = dendrocat.MasterCatalog(n1, n2, n3, catalog=irs2)
accepted = mc.catalog[mc.catalog['rejected']==0]

for i in range(len(accepted)):
    nus, fluxes, errs = mc.plotsedgrid(accepted[i], alphas=[1, 2, 3, 4], path='/users/bmcclell/nrao/documentation/SEDS_IRS2/')
    print('   nus', nus)
    print('fluxes', fluxes)
    print('errors', errs)
