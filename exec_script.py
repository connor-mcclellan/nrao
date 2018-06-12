from dendrocat import *

infilename = 'w51e2_sci.spw0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19.mfs.I.manual.image.tt0.pbcor.fits'
outfilename = 'cores.reg'

contour(infilename, outfilename, min_value=0.001, min_delta=0.002)
