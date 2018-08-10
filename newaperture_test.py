import dendrocat
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt

rs1 = dendrocat.RadioSource(fits.open('/users/bmcclell/Data/e1e2/W51e2_cont_briggsSC_tclean.image.fits.gz'))
rs1.to_catalog()
rs1.autoreject()

circ = dendrocat.aperture.Circle([40,40], 20*u.pix)
rs1.plot_grid(source_aperture=circ)
plt.show()
