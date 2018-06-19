# Source Detection, Rejection, and Flux Measurement Procedure

### Source Detection
 - Locate image file
 - Set file path to image, range of minimum values, deltas, and npix to use for source detection in "contour.py", located at the end of the script
 - Run "contour.py"
 - If plot is enabled as a keyword argument, PDF plots will be generated in "./plot_contour/"
 - A region file will be saved to "./reg/", and a catalog will be output in "./cat/"
 - Repeat for each band
 
### Rejection
 - Set the image file and catalog file locations in "reject.py", and adjust the rejection threshold if necessary
 - Run "reject.py"
 - A grid plot will show a cutout region around every detected source. Sources in grey have been rejected
 - Inspect the sources, type a manual override list
    - The list should be comma-separated. Type an "a" in front of source ID's to manually accept, and type an "r" in front of source ID's to manually reject
    - Example: to accept sources 1 and 14, and to reject sources 3 and 6, you would type "a1, a14, r3, r6". Order doesn't matter
    - The list can be typed with the plot window still open, and after closing. It will not be finalized until you press enter
 - Close the plot, and press enter to confirm the manual overrides
 - A new catalog file will be created in "./cat/" and labeled "*_filtered.cat"
 - A new region file will be created in "./reg/" and labeled "*_filtered.reg"
 - If you added manually accepted or rejected sources, run "reject.py" again so that they're saved into the catalog and region files
 - Repeat for each band
 
### Source matching
 - Add catalog filenames for each band of observation to the catalog file list in "sourcematch.py"
 - Run "sourcematch.py"
 - Sources will be matched between the provided catalog files, and a data file containing source IDs, ellipse properties, and available dendrogram fluxes across bands will be saved to "./[sky region]/mastercat_[sky_region].dat"
    - If a source is present in one image but not any others, its ellipse will still be added to the master catalog so that it can be measured across all bands for flux comparison
    - If a source is detected in multiple images, the ellipse properties will be averaged

### Flux Comparison
 - Add image filenames for each band of observation to the image file list in "fluxcomparison.py"
 - Set the master catalog file path
 - Run "fluxcomparison.py"
 - A wide range of apertures will be used for photometry across all bands
 - Fluxes in each band for each source will be added to the master catalog
