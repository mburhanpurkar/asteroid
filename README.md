# asteroid

These programs will help perform the Gaussian method of OD for an asteroid.
It is assumed that 5 star field images were taken during 3 evenly spaced
observations. From these, the orbital elements of the asteroid will be
determined.

This code was written for the. 2015 Summer Science Program at CU Boulder. If you
are a current SSPer, please do not look at this! Remember the Honor Code! Afterward, 
feel free to try it out to see how your results compare! Also, feel free to contribute! 
Most of it was written in the wee hours of the morning on a month's worth of sleep 
deprivation, so do let me know if you spot any errors or find anything that can be improved.

### LSPR
This performs centroiding and least squares plate reduction on the star field
images.

IMAGE CORRECTIONS: The star field image is first corrected using a flat image
and the background pixels are subtracted (dark current should already be
corrected for by MaxIm DL).

CENTROIDING: A weighted mean is taken using the dimensions and estimated
centroid specified by the user. The location of the centroid pixel in the image
is found.

ERROR: A variance array is generated to estimate maximum and minimum possible
centroid values. (Note that this isn't really important as the error in the orbit
determination program is what matters. This was mostly useful for developing the
code.)

LSPR: Uses least squares plate reduction to find the RA and DEC of the asteroid
given its centroid pixel values and the pixel, RA, and DEC values for six
reference stars in the star field image. (The ESO-DSS I/II image server and PPMX
optical catalogue on DS9 can be a good resource for finding reference stars.)
There is more (not at all good) error analysis code here too.

### Gaussian OD
Gauss' method of orbit determination from on  Kepler's Second Law! Outputs a nice
orbital elements chart.
