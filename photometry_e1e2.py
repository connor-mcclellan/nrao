import dendrocat
from astropy.io import fits
import astropy.units as u
from astropy.table import Table


rs1 = dendrocat.RadioSource(fits.open('/lustre/aoc/students/bmcclell/w51/W51e2_cont_briggsSC_tclean.image.fits.gz'))
rs2 = dendrocat.RadioSource(fits.open('/users/bmcclell/Data/W51e2_band3_93GHz_adjustedRADEC.fits'))
rs3 = dendrocat.RadioSource(fits.open('/users/bmcclell/Data/W51e2w_QbandAarray_adjustedRADEC.fits'))

rs3.nu = 45.0*u.GHz
rs3.freq_id = '45.0GHz'
rs3.set_metadata()

rs1.autoreject()
rs2.autoreject()

rs1.add_sources(Table.read('/users/bmcclell/nrao/cat/special_source.dat', format='ascii'))
rs1.reject([226000, 226008, 226016, 226023, 226028, 226085, 226124, 226135, 226137])
rs1.accept([226024, 226043, 226123, '226111lowd'])

rs2.accept([93005, 93006, 93007, 93012, 92021, 93054, 93053])
rs2.reject([93009, 93010, 93016, 93022, 93024, 93030, 93028, 93031, 93059])

mc = dendrocat.MasterCatalog(rs1, rs2, rs3, catalog=dendrocat.match(rs1, rs2).catalog)
mc.photometer(dendrocat.ellipse)
mc.photometer(dendrocat.annulus)

mc.catalog['_name'][mc.catalog['_name'] == '226083'] = 'w51e2e'
mc.catalog['_name'][mc.catalog['_name'] == '93045'] = 'w51e2w'
mc.catalog['_name'][mc.catalog['_name'] == '226040'] = 'w51e8'

dendrocat.utils.save_regions(mc.catalog, '/users/bmcclell/nrao/reg/mc_regions.reg')

mc.catalog.write('/users/bmcclell/nrao/cat/45-93-226GHz_photometered_adjustedRADEC.dat', format='ascii', overwrite=True)



