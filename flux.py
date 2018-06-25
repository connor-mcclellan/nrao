from glob import glob
from func import grabfileinfo
from astropy.io import fits

def flux(region):
    # import the catalog file, get names of bands
    filename = glob('mastercat_region{}*'.format(region))[0]
    bands = filename.split('bands_')[1].split('.dat')[0].split('_')
    n = len(bands)
    
    # use imgfileinfo to grab the right image files
    print("Loading image files for region {} in bands {}".format(region, bands))
    imglist = []
    for i in range(n):
        imglist.append(fits.getdata(grabfileinfo(region, bands[i])[0]))

    # set up columns for new aperture photometry measurements
    ellipse_cols = []
    circ1_cols = []
    circ2_cols = []
    circ3_cols = []
    for i in range(n)    
        ellipse_cols.append(MaskedColumn(length=len(catalog), name='ellipse_flux_band{}'.format(bands[i], mask=True))
        circ1_cols.append(MaskedColumn(length=len(catalog), name='circ1_flux_band{}'.format(bands[i], mask=True))
        circ2_cols.append(MaskedColumn(length=len(catalog), name='circ2_flux_band{}'.format(bands[i], mask=True))
        circ3_cols.append(MaskedColumn(length=len(catalog), name='circ3_flux_band{}'.format(bands[i], mask=True))
        
    # for loop iterating over each row in the catalog

        # load in this star's ellipse properties from the catalog
        
        # create all aperture shapes
        
        # measure fluxes (or noise sum within aperture for non-detection) within each shape
        
        # add fluxes to appropriate columns in the catalog

    # save catalog





