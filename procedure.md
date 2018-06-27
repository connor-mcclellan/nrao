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
 - Run "detect.py [region] [band]" with region and band arguments specified
    - ex: `python detect.py w51e2 6`
 - If plot is enabled as a keyword argument within the code, PDF plots will be generated in "./plot_contour/"
 - A region file will be saved to "./reg/", and a catalog will be output in "./cat/"
 - Repeat for each band
 
### Rejection
 - Run "reject.py [region] [band]" with region and band arguments specified
 - A grid plot will show a cutout region around every detected source. Sources in grey have been rejected
 - Inspect the sources, and type a manual override list
    - The list should be comma-separated. Type an "a" in front of source ID's to manually accept, and type an "r" in front of source IDs to manually reject
    - Example: to accept sources 301 and 314, and to reject sources 303 and 306, you would type "a301, a314, r303, r306".
 - Press enter to confirm the manual overrides
 - A new catalog file will be created in "./cat/" with "filtered" appended to the filename
 - A new region file will be created in "./reg/" with "filtered" appended to the filename
 - Repeat for each band
 
### Source matching
 - Run "match.py [region] [band1] [band2] ..." with a region argument and however many band arguments are required
    - Ex: `python match.py w51e2 3 6`
 - Sources will be matched between the catalog files with the specified region and bands, and an astropy table containing source IDs, ellipse properties, dendrogram fluxes across bands, and SNRs across bands will be saved to "./[sky region]/mastercat_region[sky_region]_bands[bands].dat"
    - If a source is present in one image but not any others, its information will still be added to the master catalog so that it can be measured across all bands for flux comparison
    - If a source is detected in multiple images, the ellipse regions will be combined
- Adjust matching threshold or perform rejection again to prevent source overlap issues

### Flux Comparison
 - Run "flux.py [region]" with a region argument specified
 - The master catalog file will be used to find positions of sources, and a wide range of apertures will be used for photometry across all bands
 - Fluxes, RMS, etc. in each band for each source will be added to the master catalog
    - Masked SNRs will be filled in at this time
    - Sources that are flagged "rejected" will remain masked, except for any fields they initially had data for
