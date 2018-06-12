from dendrocat import *

infilename = '/lustre/aoc/students/bmcclell/w51/W51e2_cont_briggsSC_tclean.image.fits.gz'

min_values = np.linspace(0.0001, 0.001, 5)
min_deltas = min_values*2.
min_npixs = [10, 15, 20]

print('Calculating dendrograms across specified parameters...')
pb = ProgressBar(len(min_values)*len(min_npixs))

for i in range(len(min_values)):
    for j in range(len(min_npixs)):
        contour(infilename, min_value=min_values[i], min_delta=min_values[i], min_npix=min_npixs[j], verbose=False)
        pb.update()
