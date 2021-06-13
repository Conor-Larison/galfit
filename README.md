# galfit
Code for modeling galaxies using GALFIT. For most updated code, go to update 3.

Basic overview: stamps.py should be run first on raw data with the mosaic Hubble image and the regions file. The regions file read-in code is hard-coded so beware if using a different image. Pretty easy fix with guess and check. This code runs through all the Lyman-alpha emitters in the field and then cuts out a 100x100 pixel (for q2343 with f140, the scale is 1 pixel = .08") image. There is a part of the code that attempts to cut out the images on the edges of the mosaic and replace the nans with random pixel values obtained from the non-nan portion of the image. This is definitely rough, but it worked for the galaxies I was worried about, so it could be improved.

Call on terminal like $python stamps.py <raw fits file> <regions file>

Next: centroids.py This code takes in the stamps and attempts to find the centroids of the centers of the stamps. The centroid is basically the flux version of center of mass, so may be around max if single galaxy, but it could be in the center of 2 - 3 galaxies if there are multiple in the center of one stamp. The box width for how the maxima are identified is hard-coded as well as the radius around the center where the centroids are allowed to reside. These can be played with to get different results from Galfit, but the ones I have in the file have worked best so far. This program also gives a flux measurment by simply taking an aperture around the centroid and summing the pixel values. From this, the code uses the photometric zeropoint (which will be different based on image so watch out), to obtain apparent magnitude measurements. This code also measures r50 and r70, which are the distances from the centroid where 50% and 70% of the flux is contained respectively. So this uses the raw flux measurement mentioned above.
  
flux.py is then called befroe gfit.py. Flux.py generates the sigma and background files that are needed for Galfit to run.
  
gfit.py This is where the actual fitting occurs. The measurements from centroids.py are read and then stored into an input file for Galfit to run. Basically the biggest things here are the different parameters, which can be found in the code. Galfit can be fastidious with what input parameters generate a good fit, so there can be a lot of trial and error here. You can also input a PSF image (generated from Tiny Tim), but the auto-generated PSF from the Galfit software seemed to work better. In the code, I have the region that is being fit defined by the locations of the centroids in the stamp and the effective radii. I found this was useful when fitting multiple-galaxy stamps, but playkng around with the maxima instead might be useful. Keeping the input image small also seemed to help the models converge with a low chi-squared (A measure of how close your model is to the raw data). The error-reading at the end of the code is extremely rough and hard-coded and can lead to a headache if the writing funciton throws an error at the end. Maybe there is an easier way to extract the uncertainties, but I am unsure.
  
gfit2.py: The same as gfit.py but a second iteration where the minimum threshold is removed so the remaining faint galaxies can also be modeled.
 
A discussion of Galfit output can be found under that related repository.

