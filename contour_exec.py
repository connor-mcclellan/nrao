from dendrocat import *

infilename = '/lustre/aoc/students/bmcclell/w51/W51e2_cont_briggsSC_tclean.image.fits.gz'

#min_values = np.linspace(0.0002, 0.00035, 2)
min_values = np.array([0.00033])
min_deltas = min_values*1.7
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
        contour(infilename, min_value=min_values[i], min_delta=min_deltas[i], min_npix=min_npixs[j], verbose=True)
        
        print('Total task progress:')
        pb.update()
