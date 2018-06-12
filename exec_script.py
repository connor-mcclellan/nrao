from dendrocat import *

infilename = '/lustre/aoc/students/bmcclell/w51/W51e2_cont_briggsSC_tclean.image.fits.gz'

contour(infilename, min_value=0.001, min_delta=0.002)
