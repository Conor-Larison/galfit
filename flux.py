# Anna de Graaff 
# Last modified 29/4/2016
#
# - Determines whether the object closest to the centre of each stamp is detected or not (above or below 3 sigma limit)
# - Plots distributions of multiples of sigma, as well as a bar chart of detected/undetected objects per field
# - Stores information (mean absolute deviation, as well as median fluxes) for use in galfit code/masking code
# - Requires output from stamps.py
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

#-------------------------------------------------------------------------------------------------------------#

basedir = '/home/conor/Project/'		# folder containing stamps
writedir = '/home/conor/Figures/'		# output folder for figures
catdir = '/home/conor/Catalogs/'		# folder where ascii files with centroids are stored
pdir = '/home/conor/Pickle/'		# output folder for pickle files
scale = 0.08 						# pixel scale [arcsec/pixel]
band = 'f140'						# select filter
field = 'q2343'				

R = int(round(0.25/scale)) # aperture radius

#-------------------------------------------------------------------------------------------------------------#

	



def open_img(filename): 
	hdulist=pf.open(filename) 
	img=hdulist[0].data
	hdulist.close()
	return img

def exp_time(filename):
	hdulist=pf.open(filename) 
	hdr=hdulist[0].header
	hdulist.close()
	return hdr['EXPTIME']

def get_name(filename):
	list = filename.split('.')
	list1 = list[0].split('/')
	obj_name = list1[-1]
	return obj_name

def circle(radius,x0,y0):
	x = np.arange(100)
	y = x.reshape(-1,1)
	d = np.sqrt((x-x0)**2+(y-y0)**2)
	mask = d<radius
	return mask

def hist1d(data,Nmin,Nmax,width):
	binl = width*np.arange(Nmin,Nmax+1)
	bincount = np.array([])
	for i in binl:
		count = np.where((data>=i)&(data<i+width))[0].size
		bincount = np.append(bincount,count)
	return width, binl, bincount

def med_abs_dev(data):	
	mdn = np.median(data)
	dev = abs(data-mdn)
	return np.median(dev)	

gt3 = np.array([])		# detected objects
lt3 = np.array([])		# undetected objects
s_total = np.array([])		# significance (in multiples of sigma/the mad)


	
	# uncomment below to single out one field
	#if field != 'q0100':
	#	continue

datadir = basedir
fits_list = sorted(glob.glob(datadir+'*.fits'))
	
d1 = {}			# store object info (significance > 3*mad, median flux, mad)
d2 = {}			# store object info (significance < 3*mad, median flux, mad)
sigma = np.array([])	# significance (in multiples of sigma/the mean absolute deviation)
mads = np.array([]) 	# median absolute deviations

	# read in centroids, select the one closest to the centre
data = ascii.read(catdir+'q2343' + '_centroids.txt')

centroids = {}	# store centroids
dc_values = {}	# store centroid distance from centre

for i, j in enumerate(data):
	obj = data['Object'][i]
	dc = data['Distance'][i]
	centroid_x = data['Centroid_X'][i]
	centroid_y = data['Centroid_Y'][i]
	if obj not in centroids:
		centroids[obj] = (centroid_x,centroid_y)
		dc_values[obj] = dc
	else: 
		if dc < dc_values[obj]:
			centroids[obj] = (centroid_x,centroid_y)


for filename in fits_list:
		
		# uncomment below to single out one file
		#if '3292' not in filename:
		#	continue

	obj = get_name(filename)
	img = open_img(filename)
	expt = exp_time(filename)
	img = img/expt

	if img.shape != (100,100): continue

		# generate random coordinates
	rdm = np.random.randint(10,90,size=(2,100))

	fluxes = np.empty([4,0])

		# calculate flux in apertures centred around random coordinates
	for i in range(len(rdm[0])):
		mask = circle(R,rdm[0][i],rdm[1][i])
		flx1 = np.sum(img[mask])
		flx2 = np.sum(img[circle(2*R,rdm[0][i],rdm[1][i])])
		flx3 = np.sum(img[circle(3*R,rdm[0][i],rdm[1][i])])
		flx4 = np.sum(img[circle(4*R,rdm[0][i],rdm[1][i])])
		flx = np.array([[flx1],[flx2-flx1],[flx3-flx2],[flx4-flx3]])
		fluxes = np.hstack((flx,fluxes))

	if obj in centroids:
		cx = centroids[obj][0]
		cy = centroids[obj][1]	
	else:
		cx = 50
		cy = 50

		# flux in aperture of increasing radius around centroid (to be used in galfit code)
	f1 = np.sum(img[circle(R,cx,cy)])	
	f2 = np.sum(img[circle(2*R,cx,cy)])
	f3 = np.sum(img[circle(3*R,cx,cy)])
	f4 = np.sum(img[circle(4*R,cx,cy)])
		# fluxes in annuli
	f21 = f2-f1
	f32 = f3-f2
	f43 = f4-f3
	f = np.array([f1,f21,f32,f43])

	s = np.array([])	# significance
	mflux = np.array([])	# median fluxes
	mad = np.array([])	# median absolute deviations for inner aperture and annuli
	for i in range(4):
		flux = f[i]
		mflux = np.append(mflux,np.median(fluxes[i]))
		mad = np.append(mad,med_abs_dev(fluxes[i]))
		s = np.append(s,abs(flux-mflux[i])/mad[i])

	mads = np.append(mads,mad[0])
	sigma = np.append(sigma,s[0])
	if s[0] >= 3: d1[obj] = (s[0],mflux,mad)
	else: d2[obj] = (s[0],mflux,mad)

d = d1.copy()
d.update(d2)
pickle.dump(d,open(pdir+field+'_sigma_'+band+'.txt','wb'))
w,l,count = hist1d(sigma,0,10,1)
plt.figure()
plt.ylim(0,max(count)+1)
plt.bar(l,count,w)
plt.ylabel('Frequency',fontsize=20)
plt.xlabel(r'$\sigma$',fontsize=20)
plt.tight_layout()
plt.savefig(writedir+field+'_sigma_'+band+'.pdf')
plt.close()	

s_total = np.append(s_total,sigma)
gt3 = np.append(gt3,len(d1))
lt3 = np.append(lt3,len(d2))
	

pickle.dump(gt3, open(pdir+'gt3_sigma_'+band+'.txt','wb'))
pickle.dump(lt3, open(pdir+'lt3_sigma_'+band+'.txt','wb'))

	
	# merge two dictionaries and store