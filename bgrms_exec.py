from dendrocat import bgrms
from astropy.io import fits

infilename = '/lustre/aoc/students/bmcclell/w51/W51e2_cont_briggsSC_tclean.image.fits.gz'
reg = './reg/reg_dend_val0.00035_delt0.000525_pix7.5.reg'

bgrms(infilename, reg)
