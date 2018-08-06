from astropy.table import MaskedColumn, Table, vstack
import astropy.units as u
import numpy as np
import time
from copy import deepcopy
from utils import commonbeam, savereg, grabcatname, grabbands
from astropy.utils.console import ProgressBar
import argparse


def dist(x, y):
    return np.sqrt(x**2 + y**2)


def getrowindex(idx, pa, t):
    """ 
    Gets the row index of a source in a table, confirming uniqueness using the _idx and the position angle.
    """
    return int(np.intersect1d(np.where(t['_idx'] == idx), np.where(t['position_angle'] == pa)))


def combine_matches(table1, table2):
    """
    Iterate through all stars in a stack of two tables, subtracting a test star's x and y from the rest of the catalog's x_cen and y_cen columns, sorting by x_cen, and testing distances starting at the top until the criterion is met."""
    complete_colnames = set(table1.colnames+table2.colnames)
    stack = vstack([table1, table2])
    stack = stack[sorted(list(complete_colnames))]
    print('Combining matches')
    numbefore = len(stack[np.where(stack['rejected']==0)])
    pb = ProgressBar(numbefore)
    i = 0
    while True:
        if i == len(stack)-1:
            break
            
        if stack[i]['rejected'] == 1:
            i += 1
            continue
            
        teststar = stack[i]
        diff_table = vstack([stack[:i], stack[i+1:]])['_idx', 'x_cen', 'y_cen', 'position_angle']
        diff_table['x_cen'] = np.abs(diff_table['x_cen'] - teststar['x_cen'])
        diff_table['y_cen'] = np.abs(diff_table['y_cen'] - teststar['y_cen'])
        diff_table.sort('x_cen')
        
        threshold = 1e-5
        found_match = False
        
        dist_col = MaskedColumn(length=len(diff_table), name='distance', mask=True)
        for j in range(10): # speed up computation by only going through 10 closest
            dist_col[j] = dist(diff_table[j]['x_cen'], diff_table[j]['y_cen'])
            if dist_col[j] <= threshold:
                found_match = True
        diff_table.add_column(dist_col)
        diff_table.sort('distance')

        if found_match:
            match_index = getrowindex(diff_table[0]['_idx'], diff_table[0]['position_angle'], stack)
            match = deepcopy(stack[match_index])
            stack.remove_row(match_index)
            
            # Find the common bounding ellipse between the match and the test star
            new_x_cen = np.average([match['x_cen'], teststar['x_cen']])
            new_y_cen = np.average([match['y_cen'], teststar['y_cen']])
            
            # REPLACE WITH COMMON BOUNDING ELLIPSE
            new_major, new_minor, new_pa = commonbeam(match['major_fwhm']*u.deg, 
                                                    match['minor_fwhm']*u.deg, 
                                                    match['position_angle']*u.deg,
                                                    teststar['major_fwhm']*u.deg,
                                                    teststar['minor_fwhm']*u.deg,
                                                    teststar['position_angle']*u.deg)

            # Save new info in test star's place
            stack[i]['x_cen'] = new_x_cen       
            stack[i]['y_cen'] = new_y_cen
            stack[i]['major_fwhm'] = new_major.value
            stack[i]['minor_fwhm'] = new_minor.value
            stack[i]['position_angle'] = new_pa.value
            
            # Replace any masked data in the teststar row with available data from the match
            for k, masked in enumerate(stack.mask[i]):       # get masked fields
                colname = stack.colnames[k]                  # get column name of masked fields
                if masked:                                   # if masked:
                    stack[i][colname] = match[colname]       # replace with data from the matched star
        i += 1
        pb.update()
    
    for colname in stack.colnames:                              # iterate over columns
        if colname.split('_')[0] == 'detected':                 # if it's a detection column
            stack[colname].fill_value = 0                       # replace masked values with 0 (False)
    
    numafter = len(stack[np.where(stack['rejected']==0)])
    print("\n{} matches combined".format(numbefore-numafter))
    return stack

    
def make_master_cat(catfilelist):
    tablelist = []
    for filename in catfilelist:
        tablelist.append(Table(Table.read(filename, format='ascii'), masked=True))
    
    current_table = tablelist[0]
    for i in range(len(tablelist)-1):
        current_table = combine_matches(current_table, tablelist[i+1])
    
    region = catfilelist[0].split('_region')[1].split('_band')[0]
    bands = grabbands(current_table)
    bandstring = ''
    for band in bands:
        bandstring += '_{}'.format(band)
    
    # Save catalog and region files
    current_table.sort('y_cen')
    current_table.write('./cat/mastercat_region{}_bands{}.dat'.format(region, bandstring), format='ascii')
    
    index = list(set(range(len(current_table)))^set(np.where(current_table['rejected']==1)[0]))
    savereg(current_table[index], './reg/masterreg_region{}_bands{}.reg'.format(region, bandstring))
        

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Match sources between bands and produce a master catalog')
    parser.add_argument('region', metavar='region', type=str, help='name of the region as listed in "imgfileinfo.dat"')
    parser.add_argument('bands', metavar='bands', type=int, nargs='+', help='integers representing the ALMA bands of observation')
    args = parser.parse_args()
    region = str(args.region)
    bands = sorted(args.bands)
    
    print("Matching sources in region {} for bands {}".format(region, bands))
    
    filelist = []
    for i in range(len(bands)):
        filelist.append(grabcatname(region, bands[i], flag='filtered'))
    
    make_master_cat(filelist)
