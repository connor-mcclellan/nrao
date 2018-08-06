from astropy.table import Table
from astropy.io import fits
import dendrocat
import astropy.units as u 
import numpy as np
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
from collections import OrderedDict

e = Table.read('/users/bmcclell/nrao/cat/w51e_external_photometered.dat', format='ascii')


rs1 = dendrocat.RadioSource(fits.open('/users/bmcclell/Data/W51e2_cont_briggsSC_tclean.image.fits.gz'))
rs2 = dendrocat.RadioSource(fits.open('/users/bmcclell/Data/W51e2_band3_93GHz_adjustedRADEC.fits'))
rs3 = dendrocat.RadioSource(fits.open('/users/bmcclell/Data/W51e2w_QbandAarray_adjustedRADEC.fits'))

rs3.nu = 45.0*u.GHz
rs3.freq_id = '45.0GHz'
rs3.set_metadata()

mc = dendrocat.MasterCatalog(rs1, rs2, rs3, catalog=e)
accepted = mc.catalog[mc.catalog['rejected']==0]

for i in range(len(accepted)):
    nus, fluxes, errs = mc.plotsedgrid(accepted[i], alphas=[1, 2, 3], 
                                       path='/users/bmcclell/nrao/documentation/SEDS_e1e2/')
    print('   nus', nus)
    print('fluxes', fluxes)
    print('errors', errs)
