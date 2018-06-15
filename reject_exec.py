from dendrocat import reject

infilename = '/lustre/aoc/students/bmcclell/w51/W51e2_cont_briggsSC_tclean.image.fits.gz'
reg = './reg/reg_dend_val0.000325_delt0.0004875_pix7.5.reg'

reject(infilename, reg)
