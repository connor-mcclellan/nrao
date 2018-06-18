## To Do
 - [X] Write code to generate a DS9 region file from a fits file, using dendrograms
 - [X] Iterate through different min_value, min_delta, min_npix settings to find best/cleanest source extraction parameters
 - [X] Define annulus regions around sources, calculate pixel RMS within those regions
 - [X] Compare annulus RMS with peak flux in the center
 - [X] Reject bad sources, create new region file with only good ones
 - [ ] Create system for tracking sources between images, so that flux measurements can be made across bands
 - [ ] Create flux histograms for all cataloged sources

## Bugs
 - SNR for noise is higher than actual sources
     - [SOLVED] When finding the peak flux, the whole image was being searched as opposed to the cutout region. Still, the image mask should have restricted those values anyway, so this might need some more looking into.
 - Lots of overlapping dendrogram regions (maybe this is ok?)
