# Source Detectionm, Rejection, and Flux Measurement Procedure

### Source Detection
 - Locate image file
 - Set range of minimum values, deltas, and npix to use for source detection in "contour.py", located at the end of the script
 - Run "contour.py"
 - PDF plots will be generated in "./plot_contour/", and a region file containing all detected sources will be saved to "./reg/"
 - Repeat for each band
 
### Rejection
 - Set the image file and region file locations in "reject.py", and adjust the rejection threshold if necessary
 - Run "reject.py"
 - A grid plot will show a cutout region around every detected source. Sources in grey have been rejected
 - Inspect the sources, and with the grid plot still open, type a comma-separated list of source numbers to manually accept the next time the script is run
 - Close the plot, and press enter.
 - A new region file will be created in "./reg/" labeled "*_filtered.reg"
 - Run "reject.py" again for the manually accepted sources to be saved into the region file
 - Repeat for each band
 
 ### Flux Measurement
 - Set the sky region name (<sky region>)
 - Set image filenames and their associated region filenames for each band of observation in "measure.py"
 - Run "measure.py"
 - Sources will be cross-matched between the three region files, and a data file containing source IDs, fluxes in each band, etc. will be saved to "<sky_region>/measurements.dat"
 - Flux histograms and source cutout plots will be saved in "./<sky region>/", labeled by source ID.
