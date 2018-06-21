from astropy.table import Table, vstack
import astropy.units as u
import numpy as np
import time
from func import convolve


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
    
    i = 0
    while True:
        print("Iteration: {}   Number of rows in stack: {}".format(i, len(stack)))
        if i == len(stack)-1:
            break
       
        teststar = stack[i]
        diff_table = vstack([stack[:i], stack[i+1:]])['_idx', 'x_cen', 'y_cen', 'position_angle']
        diff_table['x_cen'] = np.abs(diff_table['x_cen'] - teststar['x_cen'])
        diff_table['y_cen'] = np.abs(diff_table['y_cen'] - teststar['y_cen'])
        diff_table.sort('y_cen')
        

        threshold = 1e-5
        found_match = False
        d = dist(diff_table['x_cen'][0], diff_table['y_cen'][0])
        if d <= threshold:
            found_match = True

        if found_match:
            match_index = getrowindex(diff_table[0]['_idx'], diff_table[0]['position_angle'], stack)
            match = stack[match_index]
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
            for k, truth in enumerate(stack.mask[0]):
                colname = stack.colnames[k]
                if truth:
                    stack[i][colname] = match[colname]
        i += 1
        
    return stack

    
def make_master_cat(catfilelist):
    tablelist = []
    for filename in catfilelist:
        tablelist.append(Table.read(filename, format='ascii'))
    
    current_table = tablelist[0]
    for i in range(len(tablelist)-1):
        current_table = convolve_matches(current_table, tablelist[i+1])
    return current_table
        

if __name__ == '__main__':
    band3file = './cat/cat_regionw51e2_band3_val0.00015_delt0.000255_pix7.5_filtered.dat'
    band6file = './cat/cat_regionw51e2_band6_val0.000325_delt0.0005525_pix7.5_filtered.dat'
    filelist = [band3file, band6file]
    final = make_master_cat(filelist)
    print(final)
