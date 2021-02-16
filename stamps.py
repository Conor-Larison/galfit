# Anna de Graaff 
# Last modified 29/4/2016
#
# - takes HST images and reads in reg files
# - computes pixel coordinates of objects from reg file
# - create thumbnail image (stamp) of object (100x100 pixels)
#


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as col
import astropy.io.fits as pf
import os
import _pickle as pickle
from astropy import wcs
from scipy import stats
import sys
import math as mt

def main():
	data = sys.argv[1]
	reg = sys.argv[2]
	a = open(reg, 'r')
	lines = a.readlines()
	writedir = '/home/conor/Project/'

	i = 0

	while lines[i][0] != "p":
		i += 1

	lines = lines[i:]
	ra = []
	dec = []
	names = []


	for j in range(len(lines)):
		ra.append(lines[j][6:18])
		dec.append(lines[j][19:31])
		names.append(lines[j][-8:-2])
	
	a.close()




	img = img_fits(data)
	ymax = img.shape[0]
	xmax = img.shape[1]

	
	w = hdr_cood(data)
	header = get_hdr(data)
	expt = header['EXPTIME']
	obj_pix = np.empty((len(ra), 2))

	ra = hms2deg(ra)
	dec = dms2deg(dec)

	# transform coordinates and round the pixel numbers
	for k in range(len(ra)):
		world = cood_trns(w,np.array([ra[k], dec[k]]).reshape(1,2))
		obj_pix[k] = world

	new_dir = writedir

	for f in os.listdir(new_dir):	
		os.remove(os.path.join(new_dir,f))



	# create stamps
	width = 50
	for i in range(len(obj_pix)): 	
		X = obj_pix[i][0]
		Y = obj_pix[i][1]
		


		if X-width >= 0 and X+width < xmax:
			stamp = np.empty([0,2*width])
		elif X-width >= 0:
			stamp = np.empty([0,width+(xmax-X)])
		elif X+width < xmax:
			stamp = np.empty([0,width+X])
		
		for l in range(int(Y-width),int(Y+width)):
			row = np.array([])
			if l>=0 and l< ymax:
				k = int(X-width)
				while k < int(X+width):
					if k>=0 and k<xmax:
						if np.isnan(img[l][k]) == True:
							
							if k+10 >= int(X+width):
								row = np.append(row,[0]*(int(X+width) - k))
								k = int(X+width)
							else:
								row = np.append(row,[0]*10)
								k+=10
						else:
							row = np.append(row,expt*img[l][k]) # multiply pixel values by exptime for galfit to run properly
							k+=1
				stamp = np.vstack((row,stamp))

		
		xm = int(stamp.shape[1]/2)
		ym = int(stamp.shape[0]/2)



		new_stamp = []
		

		for x in range(len(stamp[0])):
			for y in range(len(stamp[1])):
				new_stamp.append(stamp[x][y])

		std = 1.4826 * stats.median_absolute_deviation(new_stamp)

		med = np.median(new_stamp)



		# check if stamp is empty and create fits file
		if np.isnan(stamp[xm][ym]) == False and stamp[xm][ym] != 0:

			for m in range(len(stamp[0])):
				for n in range(len(stamp[1])):
					
					if stamp[m][n] == 0:
						stamp[m][n] = np.random.normal(med,std)

			stamp = stamp - np.median(stamp)	# subtracts median value from image
			stamp = np.flipud(stamp)	
			hdu = pf.PrimaryHDU(stamp,header)
			hdu.writeto(new_dir+names[i]+'_q2343_f140.fits')


# functions to open fits images/headers and transform coordinates

def hdr_cood(filename):
	hdulist = pf.open(filename)
	w = wcs.WCS(hdulist[0].header)
	hdulist.close()
	return w

def get_hdr(filename):
	hdulist = pf.open(filename)
	hdr = hdulist[0].header
	hdulist.close()
	return hdr


def cood_trns(w,wcs_crd):
	pix_cood = w.wcs_world2pix(wcs_crd,1)
	return pix_cood

def hms2deg(ras):
    foo=[ra.split(':') for ra in ras]
    return np.array([float(f[0])*15+float(f[1])*15/60+float(f[2])*15/3600 for f in foo])

def dms2deg(decs):
    foo=[dec.split(':') for dec in decs]
    return np.array([float(f[0])+float(f[1])/60+float(f[2])/3600 for f in foo])


def img_fits(filename): 
	hdulist = pf.open(filename)
	img = hdulist[0].data
	hdulist.close()
	return img

main()

# Open HST images, transform wcs coordinates to pixel coordinates and cut out object stamps