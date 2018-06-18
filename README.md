## To Do
 - [X] Write code to generate a DS9 region file from a fits file, using dendrograms
 - [X] Iterate through different min_value, min_delta, min_npix settings to find best/cleanest source extraction parameters
 - [X] Define annulus regions around sources, calculate pixel RMS within those regions
 - [X] Compare annulus RMS with peak flux in the center
 - [X] Reject bad sources, create new region file with only good ones
 - [X] Constrain source detection to only dendrogram leaves
 - [X] Use dendrogram catalog instead of region file for data handling between scripts
 - [ ] Multiply ellipse dimensions by 2.35 to convert to FWHM instead of sigma
 - [ ] Add columns to dendrogram catalog for circular aperture sum, elliptical aperture sum, dendrogram contour sum
 - [ ] Use union of detected sources between images to create source IDs, so that flux measurements can be made consistently across bands
 - [ ] Create flux histograms for all cataloged sources, across all three bands

## Bugs
 - SNR for noise is higher than actual sources
     - [SOLVED] When finding the peak flux, the whole image was being searched as opposed to the cutout region. Still, the image mask should have restricted those values anyway, so this might need some more looking into.
 - Lots of overlapping dendrogram regions (maybe this is ok?)
     - [SOLVED] Use only dendrogram leaves, not branches or trunks.
