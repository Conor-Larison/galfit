# {Anna de Graaff
# Last modified 28/4/2016}
#
# Conor Larison Last modified 10/6/2021 (python 3)
#
# - Finds objects by selecting pixels that are above a certain threshold
# - Computes local maxima, and uses these to calculate object centroids
# - Creates ascii files with object number, centroids, magnitudes, fluxes and estimates for Re
# - requires output from stamps.py
#



import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as col
import astropy.io.fits as pf
from astropy.io import ascii
import glob
import os
import operator
import _pickle as pickle
from scipy.ndimage import filters
from scipy.interpolate import interp1d
from scipy import stats

#-----------------------------------------------------------------------------------------------------------#

datadir = '/home/conor/Project/'	# object stamps (obviously different based on system)
writedir = '/home/conor/Catalogs/'	# output folder
band = 'f140'					# set filter (filter used for Hubble Camera)
field = 'q2343'					# name of quasar field
scale = 0.08					# arcsec per pixel (.08 for f160, 0.04 for f475) can be found on WFC3 information on Hubble website
width = 2						# box width (hardcoded along with r0, can be played with for better results)
r0 = 10							# search radius
r1 = int(round(0.5/scale)) 			# centroid aperture radius
r2 = int(round(0.25/scale))			# flux aperture radius (used for determining significance)
zero = 26.2608				# photometric zeropoint (27.1 for f160, 25.39 for f475)
factor = 5	               # set threshold (threshold = factor*sigma) also hardcoded

#-----------------------------------------------------------------------------------------------------------#


def open_img(filename): 
	hdulist=pf.open(filename) 
	img=hdulist[0].data  
	return img

def exp_time(filename):
	hdulist=pf.open(filename) 
	hdr=hdulist[0].header 
	return hdr['EXPTIME']

def get_name(filename):
	i = 0
	while filename[i] != '_':
		i+=1
	obj_name = filename[20:i]
	return obj_name

def circle(radius,x0,y0):
	x = np.arange(100)
	y = x.reshape(-1,1)
	d = np.sqrt((x-x0)**2+(y-y0)**2)
	mask = d<radius
	return mask

def rcircle(radius,x0,y0):
	x = np.arange(100)
	y = x.reshape(-1,1)
	d = np.sqrt((x-x0)**2+(y-y0)**2)
	mask = d>radius
	return mask


	


fits_list = sorted(glob.glob(datadir+'*.fits'))

	# determine median pixel to pixel variation to use for the threshold
stdvs = np.array([])


objs = np.array([])			# object numbers
obj_x = np.array([])        # maximum x pos
obj_y = np.array([])		# maximum y pos
obj_cx = np.array([])		# centroid x pos
obj_cy = np.array([])		# centroid y pos
obj_flx = np.array([])		# flux in aperture around centroid
obj_mag = np.array([])		# object magnitude
obj_re = np.array([])		# object Re (effective radius) estimate
obj_r70 = np.array([])      # object r70 estimate
obj_std = np.array([])      # std of image
obj_offset = np.array([])   # distance from centroids to the maxima in the image





for filename in fits_list:

	img = open_img(filename)
	expt = exp_time(filename)
	for i in range(len(img)):
		if np.isnan(img[i][0]) ==  False:
			stdvs = np.append(stdvs, np.std(img[i])/expt)

	stdv = np.median(stdvs)

		
	img = open_img(filename)
	expt = exp_time(filename)
	img = img/expt
	std = 1.4826 * stats.median_absolute_deviation(img)
	simg = filters.gaussian_filter(img, 0.5) 	# create copy of image, apply Gaussian smoothing
		# select pixels that are above threshold
	threshold = factor*stdv
	rows, cols = np.where(simg > threshold)
		
	if img.shape != (100,100):
		continue


	points = {}
		# select pixels that are above threshold
	for i in range(len(rows)):
		y = rows[i]
		x = cols[i]
		if np.sqrt((50.-x)**2+(50.-y)**2) <= r0:
			points[(x,y)] = simg[y][x]

		# search for local maxima
	maxima = {}
	m = {}
	m2 = {}
	m3 = {}
		# for each point above the threshold, look for local maximum in a box of set width
	for key, value in points.items():
		d = {}
		x0 = key[0]
		y0 = key[1]
		for x in range(-width+x0,width+x0):
			for y in range(-width+y0, width+y0):
				if (x,y) in points.keys():
					d[(x,y)] = points[(x,y)]
		maxkey = max(d.items(),key=operator.itemgetter(1))[0]
		if maxkey not in m:
			m[maxkey] = max(d.values())
		# repeat procedure multiple times to converge to 1 maximum per object (only fails for large objects)
	for key1, value1 in m.items():
		d1 = {}
		x1 = key1[0]
		y1 = key1[1]
		for x in range(-width+x1,width+x1):
			for y in range(-width+y1,width+y1):
				if (x,y) in m.keys():
					d1[(x,y)] = m[(x,y)]
		maxkey1 = max(d1.items(),key=operator.itemgetter(1))[0]
		if maxkey1 not in m2:
			m2[maxkey1] = max(d1.values())

	for key2, value2 in m2.items():
		d2 = {}
		x2 = key2[0]
		y2 = key2[1]
		for x in range(-width+x2,width+x2):
			for y in range(-width+y2,width+y2):
				if (x,y) in m2.keys():
					d2[(x,y)] = m2[(x,y)]
		maxkey2 = max(d2.items(),key=operator.itemgetter(1))[0]
		if maxkey2 not in m3:
			m3[maxkey2] = max(d2.values())
		
	for key3, value3 in m3.items():
		d3 = {}
		x3 = key3[0]
		y3 = key3[1]
		for x in range(-width+x3,width+x3):
			for y in range(-width+y3,width+y3):
				if (x,y) in m3.keys():
					d3[(x,y)] = m3[(x,y)]
		maxkey3 = max(d3.items(),key=operator.itemgetter(1))[0]
		if maxkey3 not in maxima:
			maxima[maxkey3] = max(d3.values())
		
	x = np.array([])
	y = np.array([])

	
		
	for key in maxima:
		x = np.append(x,key[0])
		y = np.append(y,key[1])

	x = x.astype(int)
	y = y.astype(int)

	if len(x) > 3: #depending on search radius and box, could be more than three maxima (won't work with Galfit)
		x = np.delete(x,range(3,len(x)))
		y = np.delete(y,range(3,len(y)))

		# if no maxima are found, append (50,50)
	if len(x) == 0:
		x = np.append(x,50)
		y = np.append(y,50)


		# Find centroids
	centroidsx = np.array([])
	centroidsy = np.array([])
	flux = np.array([])
	mag = np.array([])
	obj = np.array([])
	re = np.array([])
	r70 = np.array([])
		# compute centroid
	for i in range(len(x)):
		
		weightx = 0
		weighty = 0
		flx = 0
		img_mask = rcircle(r1,x[i],y[i])
		cimg = np.copy(img)
		for j in range(x[i]-r1, x[i]+r1):
			for k in range(y[i]-r1, y[i]+r1):
				weightx = weightx + cimg[k][j]*np.float(j)
				weighty = weighty + cimg[k][j]*np.float(k)
				flx = flx + cimg[k][j]
		
		cx = weightx/flx  #It's like center of mass but for pixel values :)
		cy = weighty/flx



		
		if np.sqrt((cx-50)**2 + (cy-50)**2) <=r0: #just to make sure they're actually close to the center
			obj_offset = np.append(obj_offset,np.sqrt((cx-x[i])**2 + (cy-y[i])**2))		
			centroidsx = np.append(centroidsx, cx)
			centroidsy = np.append(centroidsy, cy)
			obj_x = np.append(obj_x,x[i])
			obj_y = np.append(obj_y,y[i])
				# store centroid
			if i == 0: obj = np.append(obj,get_name(filename))
			else: obj = np.append(obj,get_name(filename)+'_'+str(i))
				# calculate flux in small aperture (radius r2)
			f2_mask = circle(r2,x[i],y[i])
			f2 = np.sum(img[f2_mask])
			flux = np.append(flux,f2)

				# compute magnitude using larger aperture (radius r1)
			f1_mask = circle(r1,x[i],y[i])
			f1 = np.sum(img[f1_mask])
			if f1>0: mag = np.append(mag,-2.5*np.log10(f1)+zero)
			else: mag= np.append(mag,np.nan)

				# calculate flux in apertures of increasing size
			fluxes = {}
			flx70 = {}
			f_arr = np.array([])
			f_arr70 = np.array([])
			for R in range(0,r1):  #haven't looked at this in detail, but point is to find the radius at which a certain percentage of the total flux is found
				mask = circle(R,x[i],y[i])
				f = np.sum(img[mask])
				fluxes[R] = abs(f-0.5*f1)
				flx70[R] = abs(f-0.7*f1)
				f_arr = np.append(f_arr,abs(f-0.5*f1))
				f_arr70 = np.append(f_arr70,abs(f-0.7*f1))
				# estimate re using an interpolation around the absolute minimum 
			radius, value = min(fluxes.items(),key=operator.itemgetter(1))
			radius1,value1 = min(flx70.items(),key=operator.itemgetter(1))
			xmin = radius-2
			xmin1 = radius1 - 2
			if xmin < 0: xmin = 0
			if xmin1 <0: xmin1 = 0
			xmax = radius+2
			xmax1 = radius1 +2
			if xmax >= r1: xmax = r1-1
			if xmax1 >= r1: xmax1 = r1-1
			xr = np.arange(xmin,xmax+1)
			xr1 = np.arange(xmin1,xmax1+1)
			xnew = np.arange(xmin,xmax,0.1)
			xnew1 = np.arange(xmin1,xmax1,0.1)
			fc = interp1d(xr, f_arr[xmin:xmax+1],kind='quadratic')
			fc1 = interp1d(xr1, f_arr70[xmin1:xmax1+1],kind='quadratic')
			mx = np.argmin(fc(xnew))
			mx1 = np.argmin(fc1(xnew1))
			rmin = xnew[mx]
			rmin1 = xnew1[mx1]
			re = np.append(re,rmin)
			r70 = np.append(r70,rmin1)
			obj_std = np.append(obj_std,std)


	


	centroids = zip(centroidsx, centroidsy)
	objs = np.append(objs,obj)
	obj_cx = np.append(obj_cx,centroidsx)
	obj_cy = np.append(obj_cy,centroidsy)
	obj_flx = np.append(obj_flx,flux)
	obj_mag = np.append(obj_mag,mag)
	obj_re = np.append(obj_re,re)
	obj_r70 = np.append(obj_r70,r70)




dx = obj_cx - 50.		# x dev from centre
dy = obj_cy - 50.		# y dev from centre
dc = np.sqrt(dx**2 + dy**2) 	# distance from centre




	# write all data to ascii file
data = [objs, obj_x, obj_y, obj_cx, obj_cy, obj_flx, obj_mag, obj_re, obj_r70, dx, dy, dc, obj_offset]


ascii.write(data,writedir+field+'_centroids.txt',delimiter='\t',names=['Object','X','Y','Centroid_X','Centroid_Y','Flux','Magnitude','Re_est','r70 est', 'dx','dy','Distance','Offset'], overwrite = True)

