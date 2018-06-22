from astropy.table import MaskedColumn, Table, vstack
import astropy.units as u
import numpy as np
import time
from copy import deepcopy
from func import convolve, savereg, grabcatname
from astropy.utils.console import ProgressBar


def dist(x, y):
    return np.sqrt(x**2 + y**2)


def getrowindex(idx, pa, t):
    """ 
    Gets the row index of a source in a table, confirming uniqueness using the _idx and the position angle.
    """
    return int(np.intersect1d(np.where(t['_idx'] == idx), np.where(t['position_angle'] == pa)))


def convolve_matches(table1, table2):
    """
    Iterate through all stars in a stack of two tables, subtracting a test star's x and y from the rest of the catalog's x_cen and y_cen columns, sorting by y_cen, and testing distances starting at the top until the criterion is met. If the y distance exceeds the criterion, the loop terminates and no matches are found for that test star.
    """
    complete_colnames = set(table1.colnames+table2.colnames)
    stack = vstack([table1, table2])
    stack = stack[sorted(list(complete_colnames))]
    print('Convolving matches')
    pb = ProgressBar(len(stack))
    i = 0
    while True:
        if i == len(stack)-1:
            break
       
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
            
            # Convolve the match and the test star
            new_x_cen = np.average([match['x_cen'], teststar['x_cen']])
            new_y_cen = np.average([match['y_cen'], teststar['y_cen']])
            
            new_major, new_minor, new_pa = convolve(match['major_fwhm']*u.deg, 
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
            for k, truth in enumerate(stack.mask[i]):
                colname = stack.colnames[k]
                if truth:
                    stack[i][colname] = match[colname]
        i += 1
        pb.update()
        
    return stack

    
def make_master_cat(catfilelist):
    tablelist = []
    for filename in catfilelist:
        tablelist.append(Table.read(filename, format='ascii'))
    
    current_table = tablelist[0]
    for i in range(len(tablelist)-1):
        current_table = convolve_matches(current_table, tablelist[i+1])
    
    region = catfilelist[0].split('_region')[1].split('_band')[0]
    bands = [s.split('_band')[1] for s in list(filter(lambda x: x.startswith('dend_flux_band'), current_table.colnames))] # this is a disgusting hack but it works
    bandstring = ''
    for band in bands:
        bandstring += '_{}'.format(band)
    
    # Save catalog and region files
    current_table.write('./cat/mastercat_region{}_bands{}.dat'.format(region, bandstring), format='ascii')
    savereg(current_table, './reg/masterreg_region{}_bands{}.reg'.format(region, bandstring))
        

if __name__ == '__main__':
    
    region = 'w51e2'
    bands = [3, 6]
    
    filelist = []
    for i in range(len(bands)):
        filelist.append(grabcatname(region, bands[i], flag='filtered'))
    
    make_master_cat(filelist)
