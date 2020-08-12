# Anna de Graaff 
# Last modified 29/4/2016


# Further Work by Conor Larison
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
outdir = writedir + 'Output/'
fitsdir = writedir + 'Fits/'
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
#os.remove('/home/conor/Galfit/constraints.txt')
#os.remove('/home/conor/Galfit/galfit_input.txt')
#os.remove('/home/conor/Galfit/fit_results_f140.txt')

	
datadir = basedir


fits_list = sorted(glob.glob(datadir+'*.fits'))

	# read in centroids + magnitudes
data = ascii.read(catdir+field+'_centroids.txt')
centroids = {}	# centroid position for each object component
dc_values = {}	# distance of centroid from centre
	
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

filenum = 0

print(len(fits_list))
for filename in fits_list:

		
		# uncomment below to single out one object
		#if '1327' not in filename:
		#	continue

	
		
	obj = get_name(filename)
	obj_init = obj.split('_')
	obj_init = obj_init[0]


	if obj_init not in centroids:
		continue

	#if obj not in sdict: continue


	mflux = sdict[obj][1][0]
	mad = sdict[obj][2][0]

		# don't attempt to fit if object closest to centre is undetected
	if sdict[obj][0] < smin: continue

		# select centroids closest to centre (max 3)
	sort_ctds = []
	for key in dc_sorted:
		sort_ctds.append(key)



		# determine whether components are detected and set initial parameters (centroid positions, magnitude and number of components)
	if obj_init + "_1" not in sort_ctds:
		obj1 = obj_init
		comp = 1
		cx = centroids[obj1][0]
		cy = centroids[obj1][1]
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
		cx = centroids[obj1][0]
		cy = centroids[obj1][1]
		if np.isnan(centroids[obj1][2]) == False:	
			mag = centroids[obj1][2]	
		else: mag = 28	
		flx1 = centroids[obj1][3]
		cx2 = centroids[obj2][0]
		cy2 = centroids[obj2][1]
		m2 = centroids[obj2][2]
		if np.isnan(m2) == True: m2 = 28
		flx2 = centroids[obj2][3]
		s2 = (flx2-mflux)/mad
		effre2 = centroids[obj2][4]
		cx3 = centroids[obj3][0]
		cy3 = centroids[obj3][1]
		m3 = centroids[obj3][2]
		if np.isnan(m3) == True: m3 = 28
		flx3 = centroids[obj3][3]
		s3 = (flx3-mflux)/mad
		effre3 = centroids[obj3][4]
		comp = 3
	else:
		obj1 = obj_init
		obj2 = obj_init + "_1"
		cx = centroids[obj1][0]
		cy = centroids[obj1][1]
		if np.isnan(centroids[obj1][2]) == False:	
			mag = centroids[obj1][2]	
		else: mag = 28	
		flx1 = centroids[obj1][3]
		effre = centroids[obj1][4]
		cx2 = centroids[obj2][0]
		cy2 = centroids[obj2][1]
		m2 = centroids[obj2][2]
		if np.isnan(m2) == True: m2 = 28
		flx2 = centroids[obj2][3]
		s2 = (flx2-mflux)/mad
		effre2 = centroids[obj2][4]
		comp = 2
	
	comp = 1
		# create constraint file
	cf = open(writedir+'constraints.txt','w')
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
	with open(writedir+'galfit_input.txt','w') as f:
		f.write('# IMAGE and GALFIT CONTROL PARAMETERS \n')
		f.write('A) ' + filename +'\n')
		f.write('B) ' + fitsdir + 'imgblock_'+obj+'.fits \n')
		f.write('C) none \n')
		f.write('D) '+ psfdir +'\n')					# psf file
		f.write('E) ') 						# optional psf parameters'''
		f.write('F) none \n')
		f.write('H) 1 99 1 99 \n')					# image region to fit
		f.write('G) '+ writedir+'constraints.txt \n')			# constraint file
		f.write('I) 100 100 \n')					# size of convolution box
		f.write('J) '+str(zero)+' \n')					# magnitude of photometric zeropoint 25.39 for f475/ 27.1 for f160
		f.write('K) '+str(scale)+' '+str(scale)+' \n')			# plate scale
		f.write('O) regular \n')					# display type
		f.write('P) 0 \n \n')						# options 0 = normal run


		effre += 5
			# Fits first component (using Sersic parameters)
		f.write('# Object 1 \n')
		f.write(' 0) sersic \n')
		f.write(' 1) %.3f %.3f 1 1 \n' %(cx, cy)) 			# centroid position
		f.write(' 3) %.2f 1 \n' %mag)					# magnitude
		f.write(' 4) %.2f 1 \n' %effre)						# R_e
		f.write(' 5) 1. 1 \n')						# Sersic index
		f.write(' 9) 0.8. 1 \n')						# axis ratio b/a
		f.write(' 10) 0. 1 \n')					# angle of view (0 to 90 degrees)
		f.write('C0) 0. 1 \n')						# disky(-)/boxy(+)
		f.write(' Z) 0 \n \n')						# output option, 0 = residual
			
			# Fit second component
		'''if comp != 1 and s2>=smin:
			f.write('# Object 2 \n')
			f.write(' 0) sersic \n')
			f.write(' 1) %.3f %.3f 1 1 \n' %(cx2, cy2)) 		# centroid position
			f.write(' 3) %.2f 1 \n' %m2)				# magnitude
			f.write(' 4) %.2f 1 \n' %effre2)					# R_e
			f.write(' 5) 1. 1 \n')					# Sersic index
			f.write(' 9) 0.5 1 \n')					# axis ratio b/a
			f.write('C0) 0. 1 \n')					# disky(-)/boxy(+)
			f.write(' Z) 0 \n \n')					# output option, 0 = residual	
			
			# Fit third component
		if comp == 3 and s3>=smin:
			f.write('# Object 3 \n')
			f.write(' 0) sersic \n')
			f.write(' 1) %.3f %.3f 1 1 \n' %(cx3, cy3)) 		# centroid position
			f.write(' 3) %.2f 1 \n' %m3)				# magnitude
			f.write(' 4) %.2f 1 \n' %effre3)					# R_e
			f.write(' 5) 1. 1 \n')					# Sersic index
			f.write(' 9) 0.5 1 \n')					# axis ratio b/a
			f.write('C0) 0. 1 \n')					# disky(-)/boxy(+)
			f.write(' Z) 0 \n')					# output option, 0 = residual'''
			
	
	input_file = writedir+'galfit_input.txt'
		
		# checks how long fit.log file is in order to extract information from log file later
	if filenum == 0: l0 = 0
	else:
		log_stats = os.stat(outdir+'fit.log')
		l0 = log_stats.st_size		

		# Run GALFIT with the constructed input file
	subprocess.call(['/home/conor/Downloads/galfit', input_file],cwd=outdir)

		# open and store results
	result = outdir+'galfit.' + str(filenum+num).zfill(2)
	file_list = glob.glob(outdir+'galfit.*')
	if result not in file_list:
		pass
	else:	
		num += 1
		objects = np.append(objects,obj1)
		#if comp != 1  and s2>=smin: 	objects = np.append(objects,obj2)
		if comp != 1 and s3>=smin:	objects = np.append(objects,obj3)
		with open(result,'r') as f3:
			for line in f3:
				if 'Chi^2/' in line: 
					chi2 = np.append(chi2,float(line.split()[6][:-1]))
					if comp != 1 and s2>=smin:
						chi2 = np.append(chi2,float(line.split()[6][:-1]))
					if 5==4 and comp != 1 and s3>=smin:
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
		l1 = log_stats.st_size
			# extract uncertainties from log file
		with open(outdir+'fit.log','r') as f4:
			#f4.seek(-(l1-l0),2)
			lne = f4.read(l1-l0)
			keys = lne.split()
			i = keys.index('sersic')
			if len(keys[i+10]) != 1:
				i = i-1
				X_err = np.append(X_err,keys[i+11][1:-1])
			else:	X_err = np.append(X_err,keys[i+11][:-1])
			Y_err = np.append(Y_err,keys[i+12][:-1])
			mg_err = np.append(mg_err,keys[i+13])				
			Re_err = np.append(Re_err,keys[i+14])
			sc_err = np.append(sc_err,keys[i+15])
			ba_err = np.append(ba_err,keys[i+16])
				
			if 'sersic' in keys[i+18]:
				if len(keys[i+28]) != 1:
					i = i-1
					X_err = np.append(X_err,keys[i+29][1:-1])
				else:	X_err = np.append(X_err,keys[i+29][:-1])
				Y_err = np.append(Y_err,keys[i+30][:-1])
				mg_err = np.append(mg_err,keys[i+31])				
				Re_err = np.append(Re_err,keys[i+32])
				sc_err = np.append(sc_err,keys[i+33])
				ba_err = np.append(ba_err,keys[i+34])
			if len(keys) > i+36:
				if 'sersic' in keys[i+36]:
					if len(keys[i+46]) != 1:
						i = i-1
						X_err = np.append(X_err,keys[i+47][1:-1])
					else:	X_err = np.append(X_err,keys[i+47][:-1])
					Y_err = np.append(Y_err,keys[i+48][:-1])
					mg_err = np.append(mg_err,keys[i+49])				
					Re_err = np.append(Re_err,keys[i+50])
					sc_err = np.append(sc_err,keys[i+51])
					ba_err = np.append(ba_err,keys[i+52])

		
		# add constraints if fit is poor

	change = False
	if comp == 1:
		if  result not in file_list or Re[-1] > 10 or ba[-1]<0.2 or chi2[-1]>10e4:
			cf.write('1	q	0.2 to 1.0\n') 			# constrains axis ratio
			cf.write('1	re	0 to 5\n')			# constrains Re
			change = True				
	elif comp == 2:
		if result not in file_list or Re[-2] > 10 or ba[-2]<0.2 or chi2[-1]>10e4:
			cf.write('1	q	0.2 to 1.0\n') 			# constrains axis ratio for object 1
			cf.write('1	re	0 to 5\n')			# constrains Re for object 1
			change = True
		if result not in file_list or Re[-1] > 10 or ba[-1]<0.2 or chi2[-1]>10e4:
			cf.write('2	q	0.2 to 1.0\n') 			# constrains axis ratio for object 2
			cf.write('2	re	0 to 5\n')			# constrains Re	for object 2
			change = True
	'''elif comp == 3:
		if result not in file_list or Re[-3] > 10 or ba[-3]<0.2 or chi2[-1]>10e4:
			cf.write('1	q	0.2 to 1.0\n') 			# constrains axis ratio for object 1
			cf.write('1	re	0 to 5\n')			# constrains Re for object 1
			change = True
		if result not in file_list or Re[-2] > 10 or ba[-2]<0.2 or chi2[-1]>10e4:
			cf.write('2	q	0.2 to 1.0\n') 			# constrains axis ratio for object 2
			cf.write('2	re	0 to 5\n')			# constrains Re	for object 2
			change = True	
		if result not in file_list or Re[-1] > 10 or ba[-1]<0.2 or chi2[-1]>10e4:
			cf.write('3	q	0.2 to 1.0\n') 			# constrains axis ratio for object 3
			cf.write('3	re	0 to 5\n')			# constrains Re for object 3
			change = True'''
	cf.close()

	
		
		# Rerun GALFIT if input file has been altered
	'''if change == False: continue
	else: 

		subprocess.call(['/home/conor/Downloads/galfit', input_file],cwd=outdir)

			# open and store results
		result = outdir+'galfit.' + str(filenum+num).zfill(2)
		file_list = glob.glob(outdir+'galfit.*')
		if result not in file_list: continue
		else:	
			num += 1
			objects = objects[:-comp]
			X = X[:-comp]
			Y = Y[:-comp]
			Re = Re[:-comp]
			mg = mg[:-comp]
			sc = sc[:-comp]
			ba = ba[:-comp]
			X_err = X_err[:-comp]
			Y_err = Y_err[:-comp]
			Re_err = Re_err[:-comp]
			mg_err = mg_err[:-comp]
			sc_err = sc_err[:-comp]
			ba_err = ba_err[:-comp]
			chi2 = chi2[:-comp]
			objects = np.append(objects,obj1)
			if comp != 1 and s2>=smin: 	objects = np.append(objects,obj2)
			if comp != 1 and s3>=smin:	objects = np.append(objects,obj3)
			with open(result,'r') as f3:
				for line in f3:
					if 'Chi^2/' in line: 
						chi2 = np.append(chi2,float(line.split()[6][:-1]))
						if comp != 1 and s2>=smin:
							chi2 = np.append(chi2,float(line.split()[6][:-1]))
						if 5==4 and comp > 2 and s3>=smin:
							chi2 = np.append(chi2,float(line.split()[6][:-1]))
					if ' 1)' in line: 
						X = np.append(X,float(line.split()[1]))
						Y = np.append(Y,float(line.split()[2]))
					if ' 3)' in line: mg = np.append(mg, float(line.split()[1]))
					if ' 4)' in line: Re = np.append(Re,float(line.split()[1]))
					if ' 5)' in line: sc = np.append(sc,float(line.split()[1]))
					if ' 9)' in line: ba = np.append(ba, float(line.split()[1]))
			log_stats = os.stat(outdir+'fit.log')
			l2 = log_stats.st_size
			with open(outdir+'fit.log','r') as f4:
				#f4.seek(-(l2-l1),2)
				lne = f4.read(l2-l1)
				keys = lne.split()
				i = keys.index('sersic')
				if len(keys[i+10]) != 1:
					i = i-1
					X_err = np.append(X_err,keys[i+11][1:-1])
				else:	X_err = np.append(X_err,keys[i+11][:-1])
				Y_err = np.append(Y_err,keys[i+12][:-1])
				mg_err = np.append(mg_err,keys[i+13])				
				Re_err = np.append(Re_err,keys[i+14])
				sc_err = np.append(sc_err,keys[i+15])
				ba_err = np.append(ba_err,keys[i+16])
				
				if 'sersic' in keys[i+18]:
					if len(keys[i+28]) != 1:
						i = i-1
						X_err = np.append(X_err,keys[i+29][1:-1])
					else:	X_err = np.append(X_err,keys[i+29][:-1])
					Y_err = np.append(Y_err,keys[i+30][:-1])
					mg_err = np.append(mg_err,keys[i+31])				
					Re_err = np.append(Re_err,keys[i+32])
					sc_err = np.append(sc_err,keys[i+33])
					ba_err = np.append(ba_err,keys[i+34])
				if len(keys) > i+36:
					if 'sersic' in keys[i+36]:
						if len(keys[i+46]) != 1:
							i = i-1
							X_err = np.append(X_err,keys[i+47][1:-1])
						else:	X_err = np.append(X_err,keys[i+47][:-1])
						Y_err = np.append(Y_err,keys[i+48][:-1])
						mg_err = np.append(mg_err,keys[i+49])				
						Re_err = np.append(Re_err,keys[i+50])
						sc_err = np.append(sc_err,keys[i+51])
						ba_err = np.append(ba_err,keys[i+52])'''

		





# create ascii file with fit parameters

for i in range(len(objects)-len(X_err)):
	X_err = np.append(X_err,'0.0')
	Y_err = np.append(Y_err,'0.0')
	Re_err = np.append(Re_err,'0.0')
	mg_err = np.append(mg_err,'0.0')
	sc_err = np.append(sc_err,'0.0')
	ba_err = np.append(ba_err,'0.0')


output = [objects, X, X_err, Y, Y_err, Re, Re_err, mg, mg_err, sc, sc_err, ba, ba_err, chi2]


'''for array in output1:
	for i in range(10):
		if len(array[i]) == 8:
			continue
		array[i] = np.append(array[i],0)
		array = array.astype('float64')
		output.append(array)
		print(array)'''



ascii.write(output,writedir+'fit_results_'+band+'.txt', delimiter = '\t', names=['Object','Centroid_X','X_Err','Centroid_Y','Y_Err', 'R_e','R_Err', 'Magnitude','Mag_Err', 'Sersic', 'n_Err','b/a','ba_Err' ,'Chi^2'],overwrite = True)
#ascii.write(output,writedir+'fit_results_'+band+'.txt', delimiter = '\t', names=['Object'],overwrite = True)


