## To Do
 - [X] Write code to generate a DS9 region file from a fits file, using dendrograms
 - [X] Iterate through different min_value, min_delta, min_npix settings to find best/cleanest source extraction parameters
 - [X] Define annulus regions around sources, calculate pixel RMS within those regions
 - [X] Compare annulus RMS with peak flux in the center
 - [ ] Use 2x smoothed images to reject bad sources using peak flux comparison on PSFs
 - [ ] Create flux histograms for all cataloged sources

## Bugs
 - SNR for noise is higher than actual sources
 - Lots of overlapping dendrogram regions (maybe this is ok?)
