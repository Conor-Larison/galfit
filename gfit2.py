# Anna de Graaff 
# Last modified 29/4/2016
#
# - Selects (components) of objects and determines whether they are detected
# - Sets initial parameters for Galfit
# - Creates (multi-component) models to objects using Galfit
# - Outputs fits files and text files containing the fit parameters, as well as one ascii file with all parameters + uncertainties
# - Requires outputs from stamps.py, flux.py and centroids.py
#

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as col
import astropy.io.fits as pf
from astropy.io import ascii
from astropy.table import Table
import glob
import os
import subprocess
import operator
import _pickle as pickle
import shutil

#--------------------------------------------------------------------------------------------------------------------------------------#

basedir = '/home/conor/Project/'			# location of object stamps
writedir = '/home/conor/Galfit/'			# output folder for Galfit files/fits/figures
outdir = writedir + 'Output2/'
fitsdir = writedir + 'Fits2/'
figdir = writedir + 'Figures/'
catdir = '/home/conor/Catalogs/'			# location of ascii files with centroids etc
psfdir = '/home/conor/Downloads/result00_psf.fits'
pdir = '/home/conor/Pickle/'			# folder with pickle files created using flux.py code
scale = 0.08 							# arcsec per pixel, 0.08 for f160/ 0.04 for f475
band = 'f140'							# filter
zero = 26.2608							# photometric zeropoint, 25.39 for f475/ 27.1 for f160
smin = 3.							# threshold for detection of object/component
field = 'q2343'

#--------------------------------------------------------------------------------------------------------------------------------------#


def get_name(filename):
	list = filename.split(".")
	list1 = list[0].split('/')
	obj_name = list1[-1]
	return obj_name

objects = np.array([])	# object numbers
chi2 = np.array([])	# chi^2 values
X = np.array([])	# centroid x position
X_err = np.array([])	# uncertainty in x
Y = np.array([])	# centroid y position 
Y_err = np.array([])	# uncertainty in y
Re = np.array([])	# effective radius
Re_err = np.array([])	# uncertainty in Re
mg = np.array([])	# magnitudes
mg_err = np.array([])	# uncertainty in mg
sc = np.array([])	# sersic indices
sc_err = np.array([])	# uncertainty in sc
ba = np.array([])	# axis ratio b/a
ba_err = np.array([])	# uncertainty in b/a


# search for latest entry in fit.log file, gets skipped if no log file present
if not os.path.exists(outdir +'fit.log'): filenum = 0
else:
	with open(outdir+'fit.log','r') as f2:
		log_stats = os.stat(outdir+'fit.log')
		l0 = log_stats.st_size	
		if l0 < 1000:
			lne = f2.read() 
		else:		
			#f2.seek(-1000,2)
			lne = f2.read(1000)
		keys = lne.split()
		arr = np.array([])
		for key in keys:
			if 'galfit.' in key:
				arr = np.append(arr,int(key.split('.')[1]))
				filenum = int(np.max(arr))

num = 1	# dummy index to be used in for loops

# remove previous output fits files
flist = glob.glob(fitsdir+'*_'+band+'.fits')
for fits in flist:
	os.remove(fits)

shutil.rmtree(outdir)
os.mkdir(outdir)
f = open(outdir + 'fit.log','w')
f.close()



	
datadir = basedir


fits_list = sorted(glob.glob(datadir+'*.fits'))



	# read in centroids + magnitudes
names = ascii.read(writedir+ 'fit_results_f140.txt')
data = ascii.read(catdir+field+'_centroids.txt')
centroids = {}	# centroid position for each object component
dc_values = {}	# distance of centroid from centre
	

names_1 = []
for i, j in enumerate(names):
	name = names['Object'][i]
	names_1.append(name)

for i, j in enumerate(data):
	obj = data['Object'][i]
	dc = data['Distance'][i]
	centroid_x = data['Centroid_X'][i]
	centroid_y = data['Centroid_Y'][i]
	m = data['Magnitude'][i]
	flux = data['Flux'][i]
	effre = data['Re_est'][i]
	if obj not in centroids:
		centroids[obj] = (centroid_x,centroid_y,m,flux, effre)
		dc_values[obj] = dc


	
dc_sorted = sorted(dc_values,key=dc_values.get)
sdict = pickle.load(open(pdir+field+'_sigma_'+band+'.txt','rb'))


for i in names_1:
	for j in fits_list:
		if i in j:
			fits_list.remove(j)


hm = 0
filenum = 0

print(len(fits_list))

for filename in fits_list:

		
		# uncomment below to single out one object
		#if '1327' not in filename:
		#	continue

	
		
	first = filename.split('/')

	second = first[4]

	obj = second[0:6] + '_q2343_f140'

		
	obj_init = second[0:6]

	print(filename)
	print(obj_init)

	if obj_init not in centroids:
		continue

	#if obj not in sdict: continue


	mflux = sdict[obj][1][0]
	mad = sdict[obj][2][0]

		# don't attempt to fit if object closest to centre is undetected
	#if sdict[obj][0] < smin: continue

		# select centroids closest to centre (max 3)
	sort_ctds = []
	for key in dc_sorted:
		sort_ctds.append(key)



			# determine whether components are detected and set initial parameters (centroid positions, magnitude and number of components)
	if obj_init + "_1" not in sort_ctds:
		obj1 = obj_init
		comp = 1
		cx = centroids[obj1][0] +1
		cy = centroids[obj1][1] +1
		if np.isnan(centroids[obj1][2]) == False:	
			mag = centroids[obj1][2]	
		else: 
			mag = 28
		
		flx1 = centroids[obj1][3]
		effre = centroids[obj1][4]

	elif obj_init + "_2" in sort_ctds:
		obj1 = obj_init
		obj2 = obj_init + "_1"
		obj3 = obj_init + "_2"
		cx = centroids[obj1][0] + 1
		cy = centroids[obj1][1] +1
		if np.isnan(centroids[obj1][2]) == False:	
			mag = centroids[obj1][2]	
		else: mag = 28	
		flx1 = centroids[obj1][3]
		cx2 = centroids[obj2][0] +1
		cy2 = centroids[obj2][1] +1
		m2 = centroids[obj2][2]
		if np.isnan(m2) == True: m2 = 28
		flx2 = centroids[obj2][3]
		s2 = (flx2-mflux)/mad
		effre2 = centroids[obj2][4]
		cx3 = centroids[obj3][0] +1
		cy3 = centroids[obj3][1] + 1
		m3 = centroids[obj3][2]
		if np.isnan(m3) == True: m3 = 28
		flx3 = centroids[obj3][3]
		s3 = (flx3-mflux)/mad
		effre3 = centroids[obj3][4]
		comp = 3
	else:
		obj1 = obj_init
		obj2 = obj_init + "_1"
		cx = centroids[obj1][0] + 1
		cy = centroids[obj1][1] + 1
		if np.isnan(centroids[obj1][2]) == False:	
			mag = centroids[obj1][2]	
		else: mag = 28	
		flx1 = centroids[obj1][3]
		effre = centroids[obj1][4]
		cx2 = centroids[obj2][0] + 1
		cy2 = centroids[obj2][1] + 1
		m2 = centroids[obj2][2]
		if np.isnan(m2) == True: m2 = 28
		flx2 = centroids[obj2][3]
		s2 = (flx2-mflux)/mad
		effre2 = centroids[obj2][4]
		comp = 2



	#Chooses the bounds for the region to be fit based upon the locations of the centroids of the objects
	xlower = int(cx - 3*effre)
	xupper = int(cx + 3*effre)
	ylower = int(cy - 3*effre)
	yupper = int(cy + 3*effre)

	if comp ==2:
		xlower1 = int(cx - 3*effre)
		xupper1 = int(cx + 3*effre)
		ylower1 = int(cy - 3*effre)
		yupper1 = int(cy + 3*effre)

		xlower2 = int(cx2 - 3*effre2)
		xupper2 = int(cx2 + 3*effre2)
		ylower2 = int(cy2 - 3*effre2)
		yupper2 = int(cy2 + 3*effre2)

		xlower = min(xlower1,xlower2)
		xupper = max(xupper1,xupper2)
		ylower = min(ylower1,ylower2)
		yupper = max(yupper1,yupper2)
	elif comp ==3:
		xlower1 = int(cx - 3*effre)
		xupper1 = int(cx + 3*effre)
		ylower1 = int(cy - 3*effre)
		yupper1 = int(cy + 3*effre)

		xlower2 = int(cx2 - 3*effre2)
		xupper2 = int(cx2 + 3*effre2)
		ylower2 = int(cy2 - 3*effre2)
		yupper2 = int(cy2 + 3*effre2)

		xlower3 = int(cx3 - 3*effre3)
		xupper3 = int(cx3 + 3*effre3)
		ylower3 = int(cy3 - 3*effre3)
		yupper3 = int(cy3 + 3*effre3)

		xlower = min(xlower1,xlower2,xlower3)
		xupper = max(xupper1,xupper2,xupper3)
		ylower = min(ylower1,ylower2,ylower3)
		yupper = max(yupper1,yupper2,yupper3)


	effre = 2


		# create constraint file
	cf = open(writedir+'constraints2.txt','w')
	cf.write('# Component/    parameter   constraint	Comment\n# operation	(see below)   range\n')
	cf.write('1	x	-3 3\n')					# constrains x position for object 1 in a range of -3 to 3 pixels from the given centroid
	cf.write('1	y	-3 3\n')					# constrains y position for object 1 in a range of -3 to 3 pixels from the given centroid
	if comp != 1:
		cf.write('2	x	-3 3\n')				# constrains x position for object 2
		cf.write('2	y	-3 3\n')				# constrains y position for object 2
	if comp > 2:
		cf.write('3	x	-3 3\n')				# constrains x position for object 3
		cf.write('3	y	-3 3\n')				# constrains y position for object 3 

		# construct galfit input file
	with open(writedir+'galfit_input2.txt','w') as f:
		f.write('# IMAGE and GALFIT CONTROL PARAMETERS \n')
		f.write('A) ' + filename +'\n')
		f.write('B) ' + fitsdir + 'imgblock_'+obj+'.fits \n')
		f.write('C) \n')         #Sigma File
		f.write('D) \n')					# psf file
		f.write('E) \n') 						# optional psf parameters
		f.write('F) none \n')
		f.write('H) '+str(xlower)+' '+str(xupper)+' '+str(ylower)+' '+str(yupper)+' \n') # image region to fit
		f.write('G) '+ writedir+'constraints2.txt \n')			# constraint file
		f.write('I) 100 100 \n')					# size of convolution box
		f.write('J) '+str(zero)+' \n')					# magnitude of photometric zeropoint
		f.write('K) '+str(scale)+' '+str(scale)+' \n')			# plate scale
		f.write('O) regular \n')					# display type
		f.write('P) 0 \n \n')						# options 0 = normal run


			# Fits first component (using Sersic parameters)
		f.write('# Object 1 \n')
		f.write(' 0) sersic \n')
		f.write(' 1) %.3f %.3f 1 1 \n' %(cx, cy)) 			# centroid position
		f.write(' 3) %.2f 1 \n' %mag)					# magnitude
		f.write(' 4) %.2f 1 \n' %effre)						# R_e
		f.write(' 5) 1. 0 \n')						# Sersic index
		f.write(' 9) .8 1 \n')						# axis ratio b/a
		f.write(' 10) 0. 1 \n')					# angle of view (0 to 90 degrees)
		#f.write('C0) 0. 1 \n')						# disky(-)/boxy(+)
		f.write(' Z) 0 \n \n')						# output option, 0 = residual
			
			# Fit second component
		if comp != 1:
			f.write('# Object 2 \n')
			f.write(' 0) sersic \n')
			f.write(' 1) %.3f %.3f 1 1 \n' %(cx2, cy2)) 		# centroid position
			f.write(' 3) %.2f 1 \n' %m2)				# magnitude
			f.write(' 4) %.2f 1 \n' %effre2)					# R_e
			f.write(' 5) 1. 0 \n')					# Sersic index
			f.write(' 9) .8 1 \n')					# axis ratio b/a
			f.write(' 10) 0. 1 \n')					# angle of view (0 to 90 degrees)
			#f.write('C0) 0. 1 \n')					# disky(-)/boxy(+)
			f.write(' Z) 0 \n \n')					# output option, 0 = residual	
			
			# Fit third component
		if comp == 3:
			f.write('# Object 3 \n')
			f.write(' 0) sersic \n')
			f.write(' 1) %.3f %.3f 1 1 \n' %(cx3, cy3)) 		# centroid position
			f.write(' 3) %.2f 1 \n' %m3)				# magnitude
			f.write(' 4) %.2f 1 \n' %effre3)					# R_e
			f.write(' 5) 1. 0 \n')					# Sersic index
			f.write(' 9) 0.8 1 \n')					# axis ratio b/a
			f.write(' 10) 0. 1 \n')					# angle of view (0 to 90 degrees)
			#f.write('C0) 0. 1 \n')					# disky(-)/boxy(+)
			f.write(' Z) 0 \n')					# output option, 0 = residual
			
	
	input_file = writedir+'galfit_input2.txt'
		
		# checks how long fit.log file is in order to extract information from log file later
	if filenum == 0: l0 = 0
	else:
		log_stats = os.stat(outdir+'fit.log')
		l0 = log_stats.st_size + 15		

		# Run GALFIT with the constructed input file
	subprocess.call(['/home/conor/Downloads/galfit', input_file],cwd=outdir)

		# open and store results
	result = outdir+'galfit.' + str(filenum+num).zfill(2)


	file_list = glob.glob(outdir+'galfit.*')
	if result not in file_list:
		pass
	else:
		num += 1
		with open(result,'r') as f3:
			for line in f3:
				if 'Component number: 1' in line:
					objects = np.append(objects,obj1)
				if 'Component number: 2' in line:
					chi2 = np.append(chi2,0)
					objects = np.append(objects,obj2)
				if 'Component number: 3' in line:
					chi2 = np.append(chi2,0)
					objects = np.append(objects,obj3)
				if 'Chi^2/' in line: 
					chi2 = np.append(chi2,float(line.split()[6][:-1]))
				if ' 1)' in line:
					X = np.append(X,float(line.split()[1]))
					Y = np.append(Y,float(line.split()[2]))
				if ' 3)' in line: mg = np.append(mg, float(line.split()[1]))
				if ' 4)' in line: Re = np.append(Re,float(line.split()[1]))
				if ' 5)' in line: sc = np.append(sc,float(line.split()[1]))
				if ' 9)' in line: ba = np.append(ba, float(line.split()[1]))

		
			# calculate new length of log file
		log_stats = os.stat(outdir+'fit.log')
		l1 = log_stats.st_size + 15
			# extract uncertainties from log file
		with open(outdir+'fit.log','r') as f4:
			#f4.seek(-(l1-l0),2)
			
			lne = f4.read()
			keys = lne.split()
			indices = [i for i, x in enumerate(keys) if x == "sersic"]

			i = indices[-comp]


			X_err = np.append(X_err,keys[i+11][:-1])
			Y_err = np.append(Y_err,keys[i+12][:-1])
			mg_err = np.append(mg_err,keys[i+13])				
			Re_err = np.append(Re_err,keys[i+14])
			#sc_err = np.append(sc_err,keys[i+15])
			ba_err = np.append(ba_err,keys[i+16])
							
			if 'sersic' in keys[i+18]:
				X_err = np.append(X_err,keys[i+29][:-1])
				Y_err = np.append(Y_err,keys[i+30][:-1])
				mg_err = np.append(mg_err,keys[i+31])				
				Re_err = np.append(Re_err,keys[i+32])
				#sc_err = np.append(sc_err,keys[i+33])
				ba_err = np.append(ba_err,keys[i+34])
			if len(keys) > i+36:
				if 'sersic' in keys[i+36]:
					X_err = np.append(X_err,keys[i+47][:-1])
					Y_err = np.append(Y_err,keys[i+48][:-1])
					mg_err = np.append(mg_err,keys[i+49])				
					Re_err = np.append(Re_err,keys[i+50])
					#sc_err = np.append(sc_err,keys[i+51])
					ba_err = np.append(ba_err,keys[i+52])




# create ascii file with fit parameters




output = [objects,X,X_err,Y,Y_err,mg,mg_err,Re,Re_err,sc,ba,ba_err,chi2]






ascii.write(output,writedir+'fit_results_2'+band+'.txt', delimiter = '\t', names=['Object','Centroid_X','X_Err','Centroid_Y','Y_Err', 'R_e','R_Err', 'Magnitude','Mag_Err', 'Sersic','b/a','ba_Err' ,'Chi^2'],overwrite = True)
#ascii.write(output,writedir+'fit_results_'+band+'.txt', delimiter = '\t', names=['Object'],overwrite = True)




