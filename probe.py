from astropy.table import Table
from detect import detect
from func import *

infilename, region, band, _, delta_frac, _, = grabfileinfo('w51e2', 6)

min_values = np.linspace(0.00015, 0.000325, 10)
min_deltas = min_values*delta_frac
min_npixs = [5, 7.5, 10]

print("Min values: ", min_values)
print("Min deltas: ", min_deltas)
print("Min npix: ", min_npixs, '\n')

print('Total task progress:')
pb = ProgressBar(len(min_values)*len(min_npixs))
print()

for i in range(len(min_values)):
    for j in range(len(min_npixs)):
        print("\nMin value: {:.8g}    Min delta: {:.8g}   Min npix: {}".format(min_values[i], min_deltas[i], min_npixs[j]))
        detect(infilename, region, band, min_value=min_values[i], min_delta=min_deltas[i], min_npix=min_npixs[j], verbose=True)
        
        print('Total task progress:')
        pb.update()
