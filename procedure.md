# Detection, Rejection, Source Matching, and Flux Measurement Procedure

### Dendrogram Parameter Probing
 - Set the region name and band of observation
 - Set a range of minimum values, deltas, and npix to try
 - Run "probe.py"
 - Region files and catalog files for each combination of parameters will be saved to "./reg/" and "./cat/"
 - Once you find the best parameters for source detection, update the appropriate row of "imgfileinfo.dat" with the new info
 
### Source Detection
 - Optional replacement for probing if you already know specific dendrogram parameters you want to use that haven't been calculated previously
 - Update the dendrogram parameters in "imgfileinfo.dat"
 - Set region name and band of observation in "detect.py"
 - Run "detect.py"
 - If plot is enabled as a keyword argument, PDF plots will be generated in "./plot_contour/"
 - A region file will be saved to "./reg/", and a catalog will be output in "./cat/"
 - Repeat for each band
 
### Rejection
 - Set region name and band of observation in "reject.py"
 - Run "reject.py"
 - A grid plot will show a cutout region around every detected source. Sources in grey have been rejected
 - Inspect the sources, and type a manual override list
    - The list should be comma-separated. Type an "a" in front of source ID's to manually accept, and type an "r" in front of source IDs to manually reject
    - Example: to accept sources 1 and 14, and to reject sources 3 and 6, you would type "a1, a14, r3, r6". Order doesn't matter
    - The list can be typed with the plot window still open, and after closing. It will not be finalized until you press enter
 - Close the plot, and press enter to confirm the manual overrides
 - A new catalog file will be created in "./cat/" with "filtered" appended to the filename
 - A new region file will be created in "./reg/" with "filtered" appended to the filename
 - If you added manually accepted or rejected sources, run "reject.py" again so that they're saved into the catalog and region files
 - Repeat for each band
 
### Source matching
 - Set region name and bands of observation in "match.py"
 - Run "sourcematch.py"
 - Sources will be matched between the catalog files with the specified region and bands, and an astropy table containing source IDs, ellipse properties, dendrogram fluxes across bands, and SNRs across bands will be saved to "./[sky region]/mastercat_region[sky_region]_bands[bands].dat"
    - If a source is present in one image but not any others, its information will still be added to the master catalog so that it can be measured across all bands for flux comparison
    - If a source is detected in multiple images, the ellipse regions will be convolved
- Adjust matching threshold or perform rejection again to prevent source overlap issues

### Flux Comparison
 - Set region name and bands of observation in "flux.py"
 - Run "flux.py"
 - The master catalog file will be used to find positions of sources, and a wide range of apertures will be used for photometry across all bands
 - Fluxes in each band for each source will be added to the master catalog
