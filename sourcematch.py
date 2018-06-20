from astropy.table import Table


def sourcematch(catfilelist):
    catlist = []
    for catfile in catfilelist:
        catlist.append(Table.read(catfile, format='ascii'))
        
    
    
