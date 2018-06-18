# Source Detection, Rejection, and Flux Measurement Procedure

### Source Detection
 - Locate image file
 - Set file path to image, range of minimum values, deltas, and npix to use for source detection in "contour.py", located at the end of the script
 - Run "contour.py"
 - PDF plots will be generated in "./plot_contour/", a region file will be saved to "./reg/", and a catalog will be output in "./cat/"
 - Repeat for each band
 
### Rejection
 - Set the image file and catalog file locations in "reject.py", and adjust the rejection threshold if necessary
 - Run "reject.py"
 - A grid plot will show a cutout region around every detected source. Sources in grey have been rejected
 - Inspect the sources, and with the grid plot still open, type a comma-separated list of source numbers to manually accept the next time the script is run
 - Close the plot, and press enter.
 - A new catalog file will be created in "./cat/" and labeled "*_filtered.cat"
 - A new region file will be created in "./reg/" and labeled "*_filtered.reg"
 - If you added manually accepted sources, run "reject.py" again so that they're saved
 - Repeat for each band
 
 ### Flux Measurement
 - Set the sky region name ([sky region])
 - Set image filenames and their associated region filenames for each band of observation in "measure.py"
 - Run "measure.py"
 - Sources will be cross-matched between the three region files, and a data file containing source IDs, fluxes in each band, etc. will be saved to "[sky region]/measurements.dat"
