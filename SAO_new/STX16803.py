###### SAO data processing pipeline for STX16803
###### Developed by G.Lim (lim9gu@gmail.com)
###### module edited on 2020.02.29
###### run /home/lim9/anaconda3/lib/python3.7/site-packages/lgpy/SAO/STX16803.py
###### 
#=================================================================
def fileset(camera = 'STX16803'):
	"""
	1. Description
	: This function classifies calibration frames (bias, dark, skyflat, domeflat), science frames making each directory. Using other calibration frames will be added soon. If no calibration frames today, copy recent master calibration frames on current directory.
	
	2. Usage
	>>> fileset()

	3. History
	2018. 12. 28 G. Lim created.
	2019. 01. 23 Finding recent master calibration frame is added by G. Lim
	2019. 01. 24 Image size cut section is added by G. Lim
				 : Classify binned data (1x1, 2x2) with their size. Use 1x1 binned images which have 33.5MB size. 2x2 binned images are moved to bin2 directory. 
	2019. 02. 07 Loading recent calibration images is added by G. Lim
	2019. 03. 11 Image name change code is added. 
				 Domeflat, skyflat is standard naming.
				 DOMEFLAT, domeflat, SKYFLAT, Skyflat is not compatible.
	2019. 06. 06 bias, dark, flat, sci images are classified as header  information 'IMAGETYP'. and if there are no calibrations, the code calls recent master images. This function had an issue that cannot call exptime and filter which don't exist. 'hdrcheck' function is added to calculate MJD and to check if object name is entered correctly by comparing CRVAL1, CRVAL2 with IMSNG galaxy catalog. 
	2020. 02. 24 Image information is integrated as table of "cat". Remove all the space in object and filename.
	"""	
	import glob
	import os, sys
	import numpy as np
	from pyraf import iraf
	from astropy.io import fits
	from astropy.table import Table

	# Image size cut
	allimage = glob.glob('*-00*.fit')
	allimage.sort()
	xbin     = []
	ybin     = []
	IMAGETYP = []
	EXPTIME  = []
	FILTER   = []
	OBJECT   = []
	# Making image info table
	for im in allimage :
		print(im)
		hdr = fits.getheader(im)
		if (hdr['OBJECT'] == 'skyflat') | (hdr['OBJECT'] == 'Domeflat') :
			iraf.hedit(im, fields='IMAGETYP', value='Flat Field', add='yes', verify='no')
			hdr['IMAGETYP'] = 'Flat Field'
		xbin.append(hdr['XBINNING'])
		ybin.append(hdr['YBINNING'])
		IMAGETYP.append(hdr['IMAGETYP'])
		EXPTIME.append(hdr['EXPTIME'])
		if hdr['IMAGETYP'] == 'Bias Frame' :
			FILTER.append('Bias')
		elif hdr['IMAGETYP'] == 'Dark Frame' :
			FILTER.append('Dark')
		elif hdr['IMAGETYP'] == 'Flat Field' :
			FILTER.append(hdr['FILTER'])
		elif hdr['IMAGETYP'] == 'Light Frame' :
			FILTER.append(hdr['FILTER'])
		OBJECT.append(hdr['OBJECT'].replace(" ",""))
		
	cat = Table({'allimage' : allimage, 'XBINNING': xbin, 'YBINNING' : ybin,'IMAGETYP' : IMAGETYP, 'EXPTIME' : EXPTIME, 'FILTER' : FILTER, 'OBJECT' : OBJECT}, names = ['allimage', 'XBINNING', 'YBINNING', 'IMAGETYP', 'EXPTIME', 'FILTER', 'OBJECT'])

	bin1 = np.where((cat['XBINNING'] == 1) & (cat['YBINNING'] == 1))[0]
	bin2 = np.where((cat['XBINNING'] == 2) & (cat['YBINNING'] == 2))[0]
	if len(bin2) !=0 :
		print('There are 2x2 binned data...')
		os.system('/usr/bin/mkdir bin2')
		os.system('/usr/bin/mv '+" ".join(cat['allimage'][bin2])+' ./bin2')	
		print('Find 2x2 binned data at bin2 folder.')
	# Name change
	os.system('/usr/bin/rename zero bias *tion*zero.fit')
	os.system('/usr/bin/rename Bias bias *tion*Bias.fit')
	os.system('/usr/bin/rename Cal cal Cal*bias.fit')
	os.system('/usr/bin/rename D dk *tion*D*.fit')
	#os.system('rename d dk *tion*?d*.fit')
	os.system('/usr/bin/rename dark dk *tion*dark*.fit')
	os.system('/usr/bin/rename Cal cal Cal*dk*.fit')
	#calibrations = glob.glob('*.fit')
	#calibrations.sort()
	imbin1 = cat[bin1]['allimage']
	bias_idx = np.where(cat['IMAGETYP'] == 'Bias Frame')[0]
	dark_idx = np.where(cat['IMAGETYP'] == 'Dark Frame')[0]

	os.system('/usr/bin/rename Skyflat skyflat Skyflat*.fit')
	os.system('/usr/bin/rename SKYFLAT skyflat SKYFLAT*.fit')
	os.system('/usr/bin/rename SkyFlat skyflat SkyFlat*.fit')
	os.system('/usr/bin/rename domeflat Domeflat domeflat*.fit')
	os.system('/usr/bin/rename DomeFlat Domeflat DomeFlat*.fit')
	os.system('/usr/bin/rename DOMEFLAT Domeflat domeflat*.fit')
	os.system('/usr/bin/rename t-000 t_000 Domeflat*.fit')
	os.system('/usr/bin/rename flat_ flat- *flat*.fit')

	flat_idx = np.where(cat['IMAGETYP'] == 'Flat Field')[0]
	skyflat_idx = np.where((cat['IMAGETYP'] == 'Flat Field') & (cat['OBJECT'] == 'skyflat'))[0]
	domeflat_idx = np.where((cat['IMAGETYP'] == 'Flat Field') & (cat['OBJECT'] == 'domeflat'))[0]
	sci_idx  = np.where(cat['IMAGETYP'] == 'Light Frame')[0]

	bias = cat['allimage'][bias_idx]
	dark = cat['allimage'][dark_idx]
	flat = cat['allimage'][flat_idx]
	skyflat = cat['allimage'][skyflat_idx]
	domeflat = cat['allimage'][domeflat_idx]
	sci = cat['allimage'][sci_idx]
	print('All images are classified.')
	
	bias.sort()
	dark.sort()
	flat.sort()
	skyflat.sort()
	domeflat.sort()
	sci.sort()
	
	# Bias classification
	if  bias == [] :
		print('No bias today or bias directory is already created!')
		masterbias_loc = '/data1/SAO_STX16803/masterbias/'
		os.system('/usr/bin/ls '+masterbias_loc)
		masterbias_list = glob.glob(masterbias_loc+'*.fits')
		masterbias_list.sort()
		recent_masterbias = masterbias_list[-1]
		print('Recent master bias of '+recent_masterbias+' is found.')
		os.system('/usr/bin/cp '+recent_masterbias+' ./')
		print(recent_masterbias+' is copied.')
	elif bias != [] :
		os.system('/usr/bin/mkdir bias')
		biasjoin = " ".join(bias)
		os.system('/usr/bin/mv '+biasjoin+' ./bias')

	sciexptime = list(set(cat[sci_idx]['EXPTIME']))
	darkexptime = list(set(cat[dark_idx]['EXPTIME']))
	print('Check today dark exposure...')
	masterdark_loc = '/data1/SAO_STX16803/masterdark/'
	for exp in sciexptime :
		if exp in darkexptime :
			print(str(int(exp)) + ' O' )
		elif exp not in darkexptime :
			print(str(int(exp))+' X' )
			
			os.system('/usr/bin/ls '+masterdark_loc)
			masterdark_list = glob.glob(masterdark_loc+'2*dark'+str(int(exp))+'.fits')
			recent_masterdark = masterdark_list[-1]
			if len([recent_masterdark]) == 1 :
				print('Recent master dark of '+recent_masterdark+' is found.')
				os.system('/usr/bin/cp '+recent_masterdark+' ./')	
				print(recent_masterdark+' is copied.')
			elif len([recent_masterdark]) == 0 :
				print('No masterdark. Stop here :(')
	if len(dark) != 0:
		os.system('/usr/bin/mkdir dark')
		darkjoin = " ".join(dark)
		os.system('/usr/bin/mv '+darkjoin+' ./dark')	
	elif len(darkexp) == 0 :
		print('No dark today. See the copied previous master-dark.')

	# Filter of science images
	scifilter  = list(set(cat[sci_idx]['FILTER']))
	flatfilter = list(set(cat[flat_idx]['FILTER']))
	scifilter.sort()
	flatfilter.sort()
	
	for band in scifilter :
		if band in flatfilter :
			print(str(band) + ' O' )
		elif band not in flatfilter :
			print(str(band) + ' X' )
			masterflat_loc = '/data1/SAO_STX16803/masterflat_'+band+'/'
			os.system('/usr/bin/ls '+masterflat_loc)
			masterflat_list = glob.glob(masterflat_loc+'2*_n'+band+'flat.*.fits')
			recent_masterflat = masterflat_list[-1]
			if len([recent_masterflat]) == 1:
				print('Recent master flat of '+recent_masterflat+' is found.')
				os.system('/usr/bin/cp '+recent_masterflat+' ./')
			elif len([recent_masterflat]) == 0:
				print('No masterflat. Stop here :(')
	if len(skyflat) != 0:
		os.system('/usr/bin/mkdir skyflat')
		skyflatjoin = " ".join(skyflat)
		os.system('/usr/bin/mv '+skyflatjoin+' ./skyflat')	
	elif len(skyflat) == 0 :
		print('No skyflat today. See the copied previous master-flat.')

	if len(domeflat) != 0:
		os.system('/usr/bin/mkdir domeflat')
		domeflatjoin = " ".join(domeflat)
		os.system('/usr/bin/mv '+domeflatjoin+' ./domeflat')
	elif len(domeflat) == 0 :
		print('No domeflat today. See the copied previous master-flat.')			
	return bias, dark, flat, sci
#=================================================================
def biascom(camera='STX16803'):
	"""
	1. Description 
	: This function makes a master bias image of SAO 1-m using Pyraf. Put bias images to 201XXXXX/bias/ directory. Run this code on 201XXXXX directory, then pyraf chdir task enter bias directory and makes process. Output image is zero.fits and it will be copied on upper directory. Due to iraf.chdir task, you should reset python when this code is finished in order to make confusion of current directory between iraf and python! 
	
	2. Usage 
	: Start on 2018XXXX directory. Make bias directory which contains each bias frame. Naming of each bias images should be cal*bias.fit. Then just use SAO_biascom()
	
	>>> SAO_biascom()

	3. History
	2018.03	Created by G.Lim.
	2018.12.17 Edited by G.Lim. Define SAO_biascom function. 
	2019.02.07 Assign archive of masterbias in each date by G. Lim
	"""
	import glob
	import os, sys
	import numpy as np
	from pyraf import iraf
	from astropy.io import fits
	curdir = os.getcwd()
	yy   = curdir.split('/')[-1].split('-')[0]
	mm   = curdir.split('/')[-1].split('-')[1]
	dd   = curdir.split('/')[-1].split('-')[2]
	curdate = yy+mm+dd
	iraf.noao()
	iraf.imred()
	iraf.ccdred()
	iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
	iraf.chdir('./bias')
	input_name = 'bias.list'
	output_name = curdate+'_zero.fits'
	#os.system('ls cal*bias.fit > '+input_name)
	calibrations = glob.glob('cal*.fit')
	f	= open(input_name,'w+')
	for i in range(len(calibrations)):
		hdr	  = fits.getheader(calibrations[i])
		IMAGETYP = hdr['IMAGETYP']
		if IMAGETYP == 'Bias Frame' :	
			f.write(calibrations[i]+'\n')
	f.close()
	print('Zerocombine is running...')
	iraf.zerocombine(input='@'+input_name, output=output_name, combine='median', reject='minmax', process='no', scale='none', ccdtype='' )
	print('Output master '+output_name+' is created.')
	os.system('/usr/bin/cp '+output_name+' ../')
	os.system('mkdir /data1/SAO_'+camera+'/masterbias')
	os.system('/usr/bin/cp '+output_name+' /data1/SAO_'+camera+'/masterbias/')
	iraf.chdir('../')
	iraf.dir('.')
#=================================================================
def darkcom(camera='STX16803') :
	"""
	1. Description 
	: This function makes a master dark image of SAO 1-m using Pyraf. Put dark images to 201XXXXX/dark/ directory. Run this code on 201XXXXX directory, then pyraf chdir task enter dark directory and makes process. Output image is darkXXX.fits (XXX is exposure time. This function will classify each exposure of dark frames!) and it will be copied on upper directory. Due to iraf.chdir task, you should reset python when this code is finished in order to make confusion of current directory between iraf and python! 

	2. Usage
	: Start on 2018XXXX directory. Make dark directory which contains each dark frame. Naming of each dark image should be cal*dk*.fit. Then just use SAO_darkcom() 

	3. History
	2018.03	Created by G.Lim.
	2018.12.20 Edited by G.Lim. Define SAO_darkcom function.
	2019.02.07 Assign archive of masterdark in each date by G. Lim
	"""
	import glob
	import os, sys
	import numpy as np	
	from pyraf import iraf
	from astropy.io import fits
	curdir = os.getcwd()
	yy   = curdir.split('/')[-1].split('-')[0]
	mm   = curdir.split('/')[-1].split('-')[1]
	dd   = curdir.split('/')[-1].split('-')[2]
	curdate = yy+mm+dd
	iraf.noao()
	iraf.imred()
	iraf.ccdred()
	iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
	iraf.chdir('./dark')

	if camera == 'STX16803' :
		dark = glob.glob('cal*dk*.fit')
		dark.sort()
		allexptime = []
		for i in range(len(dark)) :
			hdr = fits.getheader(dark[i])
			allexptime.append(hdr['exptime'])
		expset = set(allexptime)
		exptime = list(sorted(expset))
		i=0
		for i in range(len(exptime)) :
			print('Find images with exptime of '+str(exptime[i]))
			imlist = []
			for j in range(len(dark)) :
				hdr = fits.getheader(dark[j])
				if hdr['exptime'] == exptime[i] :
					imlist.append(dark[j])
				else :
					pass
			print(imlist)
			output_name = curdate+'_dark'+str(int(exptime[i]))+'.fits'
			input_name = output_name[:-5]+'.list'
			f=open(input_name,'w+')
			for k in range(len(imlist)) : 
				f.write(imlist[k]+'\n')
			f.close()
			print('Darkcombine is running...')
			iraf.imstat('@'+input_name)
			iraf.darkcombine(input='@'+input_name, output=output_name, combine='median', reject='minmax', process='no', scale='none', ccdtype='' )
	elif camera == 'KL4040' :
		print('zero subtraction...')
		os.system('ls cal*dk*.fit > dark.list')
		iraf.imarith(operand1='@dark.list', op='-', operand2='../*_zero.fits', result='z@dark.list')
		zdark = glob.glob('z*dk*.fit')
		allexptime = []
		for i in range(len(zdark)) :
			hdr = fits.getheader(zdark[i])
			allexptime.append(hdr['exptime'])
		expset = set(allexptime)
		exptime = list(sorted(expset))
		i=0		
		for i in range(len(exptime)) :
			print('Find images with exptime of '+str(exptime[i]))
			imlist = []
			for j in range(len(zdark)) :
				hdr = fits.getheader(zdark[j])
				if hdr['exptime'] == exptime[i] :
					imlist.append(zdark[j])
				else :
					pass
			print(imlist)
			output_name = curdate+'_dark'+str(int(exptime[i]))+'.fits'
			input_name = output_name[:-5]+'.list'
			f=open(input_name,'w+')
			for k in range(len(imlist)) : 
				f.write(imlist[k]+'\n')
			f.close()
			print('Darkcombine is running...')
			iraf.imstat('@'+input_name)
			iraf.darkcombine(input='@'+input_name, output=output_name, combine='median', reject='minmax', process='no', scale='none', ccdtype='' )
	os.system('/usr/bin/cp '+output_name+' ../')
	os.system('mkdir /data1/SAO_'+camera+'/masterdark')
	os.system('/usr/bin/cp '+output_name+' /data1/SAO_'+camera+'/masterdark/')
	os.system('/usr/bin/rm d*.list')
	iraf.chdir('../')
	iraf.dir('.')
	print('Output master '+output_name+' is created.')
#=================================================================
def flatcom(flattype='sky', process='bias', camera='STX16803') :
	"""
	1. Description 
	: This function makes master-normalised images of SAO 1-m using Pyraf. Put flat images to 201XXXXX/flat/ directory. Run this code on 201XXXXX directory, then pyraf chdir task enter skyflat, domeflat directory and makes process. Keyword 'process' will decide if bias or dark subtraction is required or not. If process = bias, bias frame will be subtracted. If process = dark, dark frame will be subtracted (This function is for reducing domeflats which have long exposure time.). And then flatcombine and normalizing will be performed. Output image is nflatX.YYY.fits (X is filter and YYY is sky or dome. This function will classify each exposure of frames!) and they will be copied on upper directory. Due to iraf.chdir task, you should reset python when this code is finished in order to make confusion of current directory between iraf and python! 

	2. Usage
	: Start on 2018XXXX directory. Make skyflat or domeflat directory which contains each flat frame. Naming of each flat image should be *flat*.fit. And domeflat naming is Domeflat*.fit. Then just use SAO_flatcom(). 
	>>> SAO_flatcom('sky') --> Use skyflat
	>>> SAO_flatcom('dome') --> Use domeflat
	* default configuration is sky using bias processing.

	3. History
	2018.03	Created by G. Lim.
	2018.12.20 Edited by G. Lim. Define SAO_flatcom function.
	2018.12.28 Edited by G. Lim. Add bias or dark keyword. Join function is used when performing imarith, combine tasks.
	2019.02.07 Assign archive of masterflat in each date by G. Lim
	"""
	import glob
	import os, sys
	import itertools
	import numpy as np
	from pyraf import iraf
	from astropy.io import fits
	from astropy.io import ascii
	iraf.noao()
	iraf.imred()
	iraf.ccdred()
	iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
	curdir = os.getcwd()
	yy   = curdir.split('/')[-1].split('-')[0]
	mm   = curdir.split('/')[-1].split('-')[1]
	dd   = curdir.split('/')[-1].split('-')[2]
	curdate = yy+mm+dd
	#flattype = sys.argv[1] # dome or sky
	#if flattype == 'dome_Hal' :
	#	iraf.chdir('./domeflat_Hal')
	#elif flattype == 'dome_LED' :
	#	iraf.chdir('./domeflat_LED')
	if flattype == 'dome' :
		iraf.chdir('./domeflat')
	elif flattype == 'sky' :
		iraf.chdir('./skyflat')
	flat = glob.glob('*flat*.fit')
	flat.sort()
	if process == 'bias' : # Bias subtraction : mainly skyflat
		input_name = 'flat.list'
		k=0
		f=open(input_name,'w+')
		for k in range(len(flat)) : 
			f.write(flat[k]+'\n')
		f.close()
		#flatjoin = ",".join(flat)
		#outflatjoin = ",z".join(flat)
		#outflatjoin = 'z'+outflatjoin
		print('zero subtraction with '+flattype+'flat images...')
		iraf.imarith(operand1='@'+input_name, op='-', operand2='../*_zero.fits', result='z@'+input_name)	
		#iraf.imarith(operand1=flatjoin, op='-', operand2='../zero.fits', result=outflatjoin)	
	elif process == 'dark' : # Dark subtraction : mainly domeflat
		allexptime = []
		i=0
		for i in range(len(flat)) :
			hdr = fits.getheader(flat[i])
			allexptime.append(hdr['exptime'])
		expset = set(allexptime)
		exptime = list(sorted(expset))
		i=0
		for i in range(len(exptime)) :
			print('Find images with exptime of '+str(int(exptime[i])))
			imlist = []
			j=0
			for j in range(len(flat)) :
				hdr = fits.getheader(flat[j])
				if hdr['exptime'] == exptime[i] :
					imlist.append(flat[j])
				else :
					pass
			print(imlist)
			imlist.sort()
			input_name = 'flat.list'
			k=0
			f=open(input_name,'w+')
			for k in range(len(imlist)) : 
				f.write(imlist[k]+'\n')
			f.close()		
			# use join
			#darkjoin = ",".join(imlist)
			#outdarkjoin = ",f".join(imlist)
			#outdarkjoin = 'f'+outdarkjoin
			iraf.imarith(operand1='@'+input_name, op='-', operand2='../*_dark'+str(int(exptime[i]))+'.fits', result='d@'+input_name)
			#iraf.imarith(operand1=darkjoin, op='-', operand2='../dark'+str(int(exptime[i]))+'.fits', result=outdarkjoin)
			print('Dark subtracted flat images are created.')
	# Flat combine
	if process == 'bias' :
		calflat = glob.glob('z*.fit')
	elif process == 'dark' :
		calflat = glob.glob('d*.fit')
	allfilter = []
	i=0
	for i in range(len(calflat)) :
		hdr = fits.getheader(calflat[i])
		allfilter.append(hdr['filter'])
	filterset = set(allfilter)
	infilter = list(sorted(filterset))
	i=0
	for i in range(len(infilter)) :
		print('Find images with filter of '+str(infilter[i]))
		imlist = []
		for j in range(len(calflat)) :
			hdr = fits.getheader(calflat[j])
			if hdr['filter'] == infilter[i] :
				imlist.append(calflat[j])
			else :
				pass
		print(imlist)
		imlist.sort()
		input_name = str(infilter[i])+'flat.list'
		k=0
		f=open(input_name,'w+')
		for k in range(len(imlist)) : 
			f.write(imlist[k]+'\n')
		f.close()
		#flatjoin = ",".join(imlist)
		output_name = input_name[:-5]+'.fits'
		iraf.flatcombine(input='@'+input_name, output=output_name, combine='average', reject='crreject', process='no', scale='mode', ccdtype='', lsigma='3.', hsigma='3.' )
		#iraf.flatcombine(input=flatjoin, output=output_name, combine='median', reject='minmax', process='no', scale='mode', ccdtype='')
		print(output_name+' is created. Normalizing...')
		data, newhdr = fits.getdata(output_name, header=True)
		x = np.mean(data)
		nimage = data/x
		newflat_name = curdate+'_n'+str(infilter[i])+'flat.'+flattype+'.fits'
		fits.writeto(newflat_name, nimage, header=newhdr, overwrite=True)
		os.system('/usr/bin/cp '+newflat_name+' ../')
		os.system('mkdir '+'/data1/SAO_'+camera+'/masterflat_'+infilter[i]+'/')
		os.system('/usr/bin/cp '+newflat_name+' /data1/SAO_'+camera+'/masterflat_'+infilter[i]+'/')
	print('Normalised master flats are created.')
	iraf.imstat(images='*n?flat.'+flattype+'.fits')
	os.system('/usr/bin/rm *.list ?flat.fits')
	iraf.chdir('../')
	iraf.dir('./')
#=================================================================
def objpre(sci_list) :
	"""
	1. Description
	: This function applies master calibration images to science frames, including dark subtraction, flat fielding. STX16803 CCD has an issue of calibration frame when image type is set to bias and dark. So the CCD takes calibration images with H alpha light frame with exposure 0s as bias and >0s as dark frame.
 
	2. Usage
	: objpre('sky')

	3. History
	2018.03	Created by G. Lim
	2019.02.07 change name of master calibration frames in each date by G. Lim
	"""
	import glob
	import os, sys
	import itertools
	import numpy as np
	from pyraf import iraf
	from astropy.io import fits
	from astropy.io import ascii
	curdir = os.getcwd()
	curdate = curdir.split('/')[-1]
	iraf.noao()
	iraf.imred()
	iraf.ccdred()
	iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
	#obj_list = '*-00*.fit'
	#obj = glob.glob(obj_list)
	obj = glob.glob(sci_list)
	# dark subtraction
	allexptime = []
	i=0
	for i in range(len(obj)) :
		hdr = fits.getheader(obj[i])
		allexptime.append(hdr['exptime'])
	expset = set(allexptime)
	exptime = list(sorted(expset))
	i, j, k = 0, 0, 0
	for i in range(len(exptime)) :
		print('Find images with exptime of '+str(exptime[i]))
		imlist = []
		for j in range(len(obj)) :
			hdr = fits.getheader(obj[j])
			if hdr['exptime'] == exptime[i] :
				imlist.append(obj[j])
			else :
				pass
		print(imlist)
		imlist.sort()
		print('Creating object list for dark subtraction...')
		f = open("obj"+str(int(exptime[i]))+".list", 'w+')
		for im in range(len(imlist)) :
			f.write(imlist[im]+"\n")
		f.close()
		#imjoin = ",".join(imlist)
		#outimjoin = ",d".join(imlist)
		#outimjoin = 'd'+outimjoin
		print('dark subtraction with dark'+str(int(exptime[i]))+'.fits')
		input_dark = glob.glob('20*dark'+str(int(exptime[i]))+'.fits')[0]
		iraf.imarith(operand1 = '@obj'+str(int(exptime[i]))+'.list', op = '-', operand2 = input_dark, result = 'd@obj'+str(int(exptime[i]))+'.list')
		#iraf.imarith(operand1 = imjoin, op = '-', operand2 = './dark'+str(int(exptime[i]))+'.fits', result = outimjoin)
		print('dark subtracted object images are created.')
	# Flat fielding
	#dobj = glob.glob('d'+obj_list)
	dobj = ['d'+x for x in obj]
	dobj.sort()
	allfilter = []
	i=0
	for i in range(len(dobj)) :
		hdr = fits.getheader(dobj[i])
		allfilter.append(hdr['filter'])
	filterset = set(allfilter)
	infilter = list(sorted(filterset))
	i, j, k = 0, 0, 0
	for i in range(len(infilter)) :
		print('Find images with filter of '+str(infilter[i]))
		imlist = []
		imlist.sort()
		for j in range(len(dobj)) :
			hdr = fits.getheader(dobj[j])
			if hdr['filter'] == infilter[i] :
				imlist.append(dobj[j])
			else :
				pass
		print(imlist)

		g = open("obj"+infilter[i]+".list", 'w+')
		for im in range(len(imlist)) :
			g.write(imlist[im]+"\n")
		g.close()
		#dimjoin = ",".join(imlist)
		#doutimjoin = ",f".join(imlist)
		#doutimjoin = 'f'+doutimjoin
		print('Performing flat fielding...')
		#nflats = glob.glob('2*n'+str(infilter[i])+'flat.'+flattype+'.fits')[0]
		nflats = glob.glob('2*n'+str(infilter[i])+'flat.*.fits')[0]
		flattype = nflats.split('.')[1]
		print(nflats)
		#nflats = 'n'+str(infilter[i])+'flat.'+flattype+'.fits'
		iraf.imarith(operand1='@obj'+infilter[i]+'.list', op='/', operand2=nflats, result='f@obj'+infilter[i]+'.list')
		#iraf.imarith(operand1=dimjoin, op='/', operand2=nflats, result=doutimjoin)
	print('Flat fielding is finished. Check the images.')
#=================================================================
def astrometry() :
	"""
	1. Description 
	: Solving WCS coordinates using Astrometry.net software. For better performance in especially B band images, --use-sextractor mode is added. This mode needs SExtractor configuration files. So please posit configuration files for your working directory. cpulimit 300 is also added to prevent too long processing time for bad images.  

	2. Usage
	>>> SAO_astrometry()

	3. History
	2018.03	Created by G.Lim.
	2018.12.18 Edited by G.Lim. SExtractor mode is added.
	2018.12.21 Edited by G.Lim. Define SAO_astrometry function.
	"""
	import os,sys
	import glob
	import subprocess
	import numpy as np
	addlist = glob.glob('fd*.fit')
	addlist.sort()
	sexconfig = '/data1/code/astrom.config/astrometry.net.sex'
	print('Solving WCS using Astrometry.net...')
	for n in range(len(addlist)):
		#com='/usr/bin/solve-field '+addlist[n]+' --cpulimit 300 --overwrite --use-sextractor  --sextractor-config '+sexconfig+' --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_AUTO --sort-ascending --scale-unit arcsecperpix --scale-low 0.2 --scale-high 0.4 --no-remove-lines --uniformize 0 --no-plots  --new-fits a'+addlist[n]+' --temp-dir .\n'
		# for RASA36
		#com='/usr/bin/solve-field '+addlist[n]+' --cpulimit 300 --overwrite --use-sextractor  --sextractor-config '+sexconfig+' --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_AUTO --sort-ascending --scale-unit arcsecperpix --scale-low 2.8 --scale-high 3 --no-remove-lines --uniformize 0 --no-plots  --new-fits a'+addlist[n]+' --temp-dir .\n'
		# for SAO KL4040
		com='/usr/bin/solve-field '+addlist[n]+' --cpulimit 300 --overwrite --use-sextractor  --sextractor-config '+sexconfig+' --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_AUTO --sort-ascending --scale-unit arcsecperpix --scale-low 0.29 --scale-high 0.32 --no-remove-lines --uniformize 0 --no-plots  --new-fits a'+addlist[n]+' --temp-dir .\n'
		print(com)
		print(str(n)+' th of '+str(len(addlist)))
		os.system(com) 
	orinum = subprocess.check_output('ls fd*.fit | wc -l', shell=True)
	resnum = subprocess.check_output('ls afd*.fit | wc -l', shell=True)
	print("from "+str(orinum[:-1])+" files , "+str(resnum[:-1])+" files are solved.")
	print("All done.")
	os.system('rm tmp*')
	os.system('rm *.wcs *.rdls *.corr *.xyls *.solved *.axy *.match ')
	print('Astrometry process is complete.')
#=================================================================
def hdrcheck(imlist_name):
	"""
	Check header and put galaxy name comparing to IMSNG catalog. 
	"""
	import os
	import sys
	import glob
	import astropy.units as u
	from astropy.time import Time
	from astropy.io import ascii
	from astropy.io import fits
	from astropy.coordinates import SkyCoord

	SAO_fov = 21.2 # arcmin
	all_catname = '/data1/code/alltarget.dat'
	all_cat	 = ascii.read(all_catname) 
	ra, dec = all_cat['ra'], all_cat['dec']
	radeg, decdeg = [], []
	for i in range(len(all_cat)) :
		c   = SkyCoord(str(ra[i])+' '+str(dec[i]), unit=(u.hourangle, u.deg))
		radeg.append(c.ra.deg)
		decdeg.append(c.dec.deg)
	all_cat['radeg'] = radeg
	all_cat['decdeg'] = decdeg
	coo_all = SkyCoord(radeg, decdeg, unit=(u.deg,u.deg))


	imlist = glob.glob(imlist_name)
	imlist.sort()
	for i in range(len(imlist)) :
		inim			= imlist[i]
		print(inim)
		data, hdr	   = fits.getdata(inim, header=True)
		CRVAL1		  = hdr['CRVAL1']
		CRVAL2		  = hdr['CRVAL2']
		camera = inim.split('-')[1]
		mjd0 = 2400000.5
		if camera == 'MAO_SNUCAM' :
			t = Time(hdr['utdate']+'T'+hdr['utstart'], format='isot', scale='utc')
			jd = t.jd
			mjd = t.mjd
			hdr['JD'] = round(jd,5)
		else: 
			jd = hdr['jd']
			mjd  = jd - mjd0
		hdr['MJD']  = round(mjd,5)
		coo_target	  = SkyCoord(CRVAL1, CRVAL2, unit=(u.deg, u.deg))
		indx, d2d, d3d  = coo_target.match_to_catalog_sky(coo_all)

		if d2d.arcmin > SAO_fov/2. :
			print('Coordinates of the image are not in IMSNG catalog. No matching. Maybe you obtained wrong field. OR Non-IMSNG target.' )
			fits.writeto(inim, data, header=hdr, overwrite=True)
			print('Only MJD is entered in image header.')
			pass
		elif d2d.arcmin < SAO_fov/2. :
				obj  = all_cat[indx]['obj']
				print('======================================')
				print(obj+ ' is matched.')
				print(str(round(d2d.arcmin[0],3))+ ' arcmin apart')
				print('======================================')
				hdr['object']   = obj
				fits.writeto(inim, data, header=hdr, overwrite=True)
				
	print('Header info inspection is finished.')
#=================================================================
def fnamechange(imlist_name,camera='SAO_STX16803') :
	"""
	1. Description 
	: Change file name of WCS solving images (afd*.fits) using naming sequence following: Calib-SAO-CCD_name-OBJECT-UTDATE-UTSTART-FILTER-EXPTIME.fits. Then the code copies changed files to IMSNGgalaxies directory.

	2. Usage
	>>> SAO_fnamechange()

	3. History
	2018.03	Created by G.Lim.
	2018.12.21 Edited by G.Lim. Define SAO_fnamechange function.
	"""
	import glob
	import os, sys
	import subprocess
	import numpy as np 
	from astropy.io import fits
	lists = glob.glob(imlist_name)
	lists.sort()
	for i in range(len(lists)):
		hdr = fits.getheader(lists[i])
		EXPTIME = int(hdr['exptime'])
		if camera == 'Kepler' :
			FILTER = 'R'
		elif camera == 'SAO_STX16803' :
			FILTER = hdr['filter']
		elif camera == 'CCA250' :
			OBJECT = hdr['object']
			FILTER = hdr['filter']
		if (OBJECT == 'M51') | (OBJECT == 'M51a') :
			print('Object name is corrected.')
			OBJECT = 'M51A'
		UTDATE = hdr['date-obs'][0:10]
		UTSTART = hdr['date-obs'][11:19]
		newimage = 'Calib-'+camera+'-'+OBJECT+'-'+str(UTDATE[0:4])+str(UTDATE[5:7])+str(UTDATE[8:10])+'-'+str(UTSTART[0:2])+str(UTSTART[3:5])+str(UTSTART[6:8])+'-'+FILTER+'-'+str(EXPTIME)+'.fits'
		os.system('/usr/bin/cp '+lists[i]+' '+newimage)
		print('Copy '+lists[i]+' to '+newimage+'.')
	#os.system('cp Cal*.fits /data3/IMSNG/IMSNGgalaxies/')
	print("Basic preprocessing of SAO 1-m is finished.")
#=================================================================
def filemove(camera='STX16803'):
	"""
	1. Description 
	: This code distributes all the calibrated images of SAO 1m data to IMSNGgalaxies/SAO directory, based on C. Choi's code.
	
	2. Usage
	: Run this code on '/data3/IMSNG/IMSNGgalaxies' location.
	>>> SAO_filemove() 

	3. History
	2018.12	Created by G.Lim 
	2018.01.24 Docstring is added by G. Lim
	"""
	import os
	import numpy as np
	import astropy.io.fits as fits
	import astropy.io.ascii as ascii
	loc = '/data3/IMSNG/IMSNGgalaxies/'
	os.system('ls Cal*.fits -l > currentdir.list')
	curdirlist = ascii.read('currentdir.list',data_start=0)
	dirpar	 = curdirlist['col2']
	name	   = curdirlist['col9']
	#dirnames   = name[np.where(dirpar >  1)]
	filenames  = name[np.where(dirpar == 1)]
	calframes  = []
	for i in filenames : 
		if (i[:9] =='Calib-SAO') & (i[-4:] == 'fits') : calframes.append(i)
	print(str(len(calframes))+' files exist.')
	for n in calframes:
		galname  = n.split('-')[2]
		print(galname)
		if galname == 'M51a' :
			galname = 'M51A'
		hdr = fits.getheader(n)
		if camera == 'STX16803' :
			infilter = hdr['filter']
		elif camera == 'Kepler' :
			infilter = 'R'
		makedir	= '/usr/bin/mkdir '+loc+galname
		os.system(makedir)
		makeobsdir = '/usr/bin/mkdir '+loc+galname+'/SAO/' # If IMSNG target
		os.system(makeobsdir)
		makecamdir = '/usr/bin/mkdir '+loc+galname+'/SAO/'+camera
		os.system(makecamdir)
		makefildir = '/usr/bin/mkdir '+loc+galname+'/SAO/'+camera+'/'+infilter+'/'
		os.system(makefildir)
		mvcommand='/usr/bin/cp '+ n +' '+loc+galname+'/SAO/'+camera+'/'+infilter+'/'
		mvsubtract = '/usr/bin/cp hd'+n+' hc'+n+' '+loc+galname+'/SAO/'+camera+'/'+infilter+'/'
		os.system(mvcommand)
		os.system(mvsubtract)
		chmodcom = '/usr/bin/chmod 777 -R '+loc+galname+'/SAO/'
		os.system(chmodcom)
#=================================================================
def makedir(camera='STX16803'):
	"""
	1. Description
	: Make directory of camera name (STX16803 or Kepler)

	2. Usage
	>>> 

	3. History
	2019.04.21 : G. Lim created.
	"""
	import glob
	import os, sys
	workdir = '/data3/IMSNG/IMSNGgalaxies'
	curdir = os.getcwd()
	if curdir != workdir :
		print('Please run this code in '+workdir)
	elif curdir == workdir :
		folder = glob.glob('*')
		folder.sort()
		for i in range(len(folder)):
			os.chdir(folder[i])
			obs = glob.glob('*')
			obs.sort()
			if 'SAO' in obs :
				print('No camera dir. Make '+camera)
				os.chdir('SAO')
				os.system('/usr/bin/mkdir '+'./'+ camera)
				os.system('/usr/bin/mv Cal*.fits '+camera)
				#os.system('rsync -a B '+camera+'/B')
				os.chdir('../../')
			else :
				print(folder[i]+' has no SAO data.')
				os.chdir('../')
				pass
		print('Done.')
#=================================================================
def run_gregister_list(imlist_name) :
	"""
	1. Description 
	: 
	ref = True : specify
	2. Usage
	>>> SAO_gregister() 

	3. History
	2019.04.17  Created by G. Lim
	"""
	import glob
	import alipy
	import os,sys
	import numpy as np
	from astropy.io import fits

	### Object classification
	# imlist_name = 'Cal*skysub.fits'
	lists	  = glob.glob(imlist_name)
	lists.sort()
	split = lists[0].split('-')
	objects	= []
	for i in range(len(lists)) :
		hdr	= fits.getheader(lists[i])
		objects.append(hdr['object'])
	objectset  = set(objects)
	objectlist = list(sorted(objectset))

	def gregister(images):
		id	 = alipy.ident.run(ref_image, images, visu=False)
		print("%20s : %20s, flux ratio %.2f" % (id[0].ukn.name, id[0].trans, id[0].medfluxratio))
		alipy.align.irafalign(id[0].ukn.filepath, id[0].uknmatchstars, id[0].refmatchstars, shape=outputshape, makepng=False)

	### Filter classification
	obj = 0
	for obj in range(len(objectlist)):
		object_name = objectlist[obj]
		if 'gregister' in split :
			image_list  = glob.glob(split[0]+'*'+object_name+'*gre*.fits')
		elif 'skysub' in split :
			image_list  = glob.glob(split[0]+'*'+object_name+'*skysub.fits')
		elif ('gregister' in split) & ('skysub' in split) :
			image_list  = glob.glob(split[0]+'*'+object_name+'*skysub*gre*.fits')
		elif ('gregister' not in split) & ('skysub' not in split) :
			image_list  = glob.glob(split[0]+'*'+object_name+'*0.fits')
		elif len(objectlist) == 1 :
			image_list  = objectlist[0]
		#image_list  = glob.glob('Cal*'+object_name+'*0.fits')
		k		   = 0
		allfilter   = []
		for k in range(len(image_list)) :
			hdr	 = fits.getheader(image_list[k])
			allfilter.append(hdr['filter'])
		filterset   = set(allfilter)
		infilter	= list(sorted(filterset))	
		band		= 0
		for band in range(len(infilter)):
			if 'gregister' in split :
				image_list_filter	 = glob.glob(split[0]+'*'+object_name+'*'+infilter[band]+'*gre*.fits')
			elif 'skysub' in split :
				image_list_filter	 = glob.glob(split[0]+'*'+object_name+'*'+infilter[band]+'*skysub.fits')
			elif ('gregister' in split) & ('skysub' in split) :
				image_list_filter	 = glob.glob(split[0]+'*'+object_name+'*'+infilter[band]+'*skysub*gre*.fits')
			elif ('gregister' not in split) & ('skysub' not in split) :
				image_list_filter	 = glob.glob(split[0]+'*'+object_name+'*'+infilter[band]+'*0.fits')
			if len(image_list_filter) > 1 : 
				images_to_align   = image_list_filter
				ref_image = image_list_filter[1]
				outputshape	   = alipy.align.shape(ref_image)
				images_to_align_1 = []
				n = 0
				for n in range(len(images_to_align)) :
					images		= images_to_align[n:n+1][0]
					print('\n',images,'\n')
					gregister([images])
					images_to_align_1.append(images)
			else :
				print(str(image_list_filter[0])+' is the only element. Pass.')
				pass	
	os.system('mv alipy_out/*.fits .')
	os.system('rm -r alipy_out')
	print('Done. \a')
#=================================================================
def gregistering(images_to_align, ref_image):
	import os, sys, glob
	import alipy
	import time
	from pyraf import iraf
	import numpy as np
	from astropy.io import fits
	import astropy.units as u
	from astropy.time import Time
	starttime = time.time()
	if ref_image == '':
		ref_image = images_to_align[0]
	identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
	for id in identifications:  # list of the same length as images_to_align.
		if id.ok == True:  # i.e., if it worked
			print("%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio))
		else:
			print("%20s : no transformation found !" % (id.ukn.name))
	outputshape = alipy.align.shape(ref_image)
	for id in identifications:
		if id.ok == True:
			params_align = dict(	filepath=id.ukn.filepath,
                            uknstarlist=id.uknmatchstars,
                            refstarlist=id.refmatchstars,
                            shape=alipy.align.shape(ref_image),
                            outdir='./',
                            makepng=False)
			alipy.align.irafalign(**params_align)
	deltime = time.time() - starttime
	print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')
#=================================================================
def run_gregister_epoch(imlist_name, sep=10.):
	import glob
	import alipy
	import os,sys
	import numpy as np
	from astropy.io import fits
	from lgpy.SameEpoch import SameEpoch	
	from lgpy.SAO import STX16803
	res = SameEpoch(imlist_name, sep=sep)
	for i in range(len(res)):
		#imlist = list(res[i])
		STX16803.gregistering([list(res[i])][0],[list(res[i])][0][0])
#==================================================================
def run_gregister(imlist_name, refim, input_list=False, large=True):
	"""
	1. Description 
	: Specify reference image into refim. images to align are imlist_name.

	2. Usage
	>>> gregister_ref('Cal*.fits', 'ref', large=True) 

	large = True : Sci's pixscale is the same with or larger than refim (sci images are processed.)
	large = False : Sci's pixscale is smaller than refim.
	(ref image will be copied and processed.)

	3. History
	2019.06.13  Created by G. Lim
	"""
	import glob
	import alipy
	import os,sys
	import numpy as np
	from astropy.io import fits

	def gregister(ref, inim):
		id	 = alipy.ident.run(ref, inim, visu=False)
		print("%20s : %20s, flux ratio %.2f" % (id[0].ukn.name, id[0].trans, id[0].medfluxratio))
		alipy.align.irafalign(id[0].ukn.filepath, id[0].uknmatchstars, id[0].refmatchstars, shape=outputshape, makepng=False)	
		return [id[0].ukn.filepath]
	if input_list == False :
		image_list_filter	  = glob.glob(imlist_name)	
	elif input_list == True :
		image_list_filter = list(imlist_name)
	image_list_filter.sort()

	if large == False :
		if len(image_list_filter) > 1 : 
			images_to_align   = image_list_filter
			print('ref : '+refim)
			outputshape	   = alipy.align.shape(refim)
			images_to_align_1 = []
			n = 0
			for n in range(len(images_to_align)) :
				images		= images_to_align[n:n+1][0]
				#images		= images_to_align
				print('\n',images,'\n')
				gregister(refim, images)
				images_to_align_1.append(images)
				outim = id.ukn.filepath[:-5]+'_gregister.fits'
				newoutim = 'g'+id.ukn.filepath[:-5]+'.fits'
				os.system('cp ./alipy_out/'+outim+' ./'+newoutim)
		else :
			print(str(image_list_filter[0])+' is the only element. Pass.')
			pass

	elif large == True :
		if len(image_list_filter) > 1 : 
			for i in range(len(image_list_filter)):
				inim =  image_list_filter[i]
				part = inim.split('-')
				part[0] = 'Ref'
				newref = '-'.join(part)
				os.system('cp '+refim+' '+newref)
				images_to_align = newref
				refim= inim
				print('ref : '+inim) 
				outputshape	   = alipy.align.shape(refim)
				images_to_align_1 = []
				images		= images_to_align
				gregister([images])
				images_to_align_1.append(images)
				outim = id.ukn.filepath[:-5]+'_gregister.fits'
				newoutim = 'g'+id.ukn.filepath[:-5]+'.fits'
				os.system('cp ./alipy_out/'+outim+' ./'+newoutim)
		else :
			print(str(image_list_filter[0])+' is the only element. Pass.')
			pass
	os.system('mv alipy_out/*.fits .')
	os.system('rm -r alipy_out')
	print('Done. \a')
#=================================================================
'''
def run_gregister_epoch(imlist_name, sep=10.):
	import glob
	import alipy
	import os,sys
	import numpy as np
	from astropy.io import fits
	from lgpy.SameEpoch import SameEpoch	
	from lgpy.SAO import STX16803
	res = SameEpoch(imlist_name, sep=sep)
	for i in range(len(res)):
		#imlist = list(res[i])
		STX16803.run_gregister(res[i], str(res[i][0]), input_list=True, large=False)
'''
#=================================================================
'''
def run_gregister_epoch(imlist_name, large=True):
	"""
	1. Description 
	: Specify reference image into refim. images to align are imlist_name.

	2. Usage
	>>> gregister_ref('Cal*.fits', 'ref', large=True) 

	large = True : Sci's pixscale is the same with or larger than refim (sci images are processed.)
	large = False : Sci's pixscale is smaller than refim.
	(ref image will be copied and processed.)

	3. History
	2019.06.13  Created by G. Lim
	"""
	import glob
	import alipy
	import os,sys
	import numpy as np
	from astropy.io import fits
	from lgpy.SameEpoch import SameEpoch
	

	def gregister(images):
		id	 = alipy.ident.run(ref_image, images, hdu=0, visu=False)
		print("%20s : %20s, flux ratio %.2f" % (id[0].ukn.name, id[0].trans, id[0].medfluxratio))
		alipy.align.irafalign(id[0].ukn.filepath, id[0].uknmatchstars, id[0].refmatchstars, shape=outputshape, makepng=False)	
		return [id[0].ukn.filepath]

	res = SameEpoch(imlist_name, sep=sep)
	for i in range(len(res)):
		image_list_filter = list(res[i])
		image_list_filter.sort()

		if large == False :
			if len(image_list_filter) > 1 : 
				images_to_align   = image_list_filter
				ref_image		  = image_list_filter[0]
				print('ref : '+ref_image)
				outputshape	      = alipy.align.shape(ref_image)
				images_to_align_1 = []
				n = 0
				for n in range(len(images_to_align)) :
					images		= images_to_align[n:n+1][0]
					print('\n',images,'\n')
					idfile = gregister([images])
					images_to_align_1.append(images)
					outim = idfile[0][:-5]+'_gregister.fits'
					newoutim = 'g'+idfile[0][:-5]+'.fits'
					os.system('cp ./alipy_out/'+outim+' ./'+newoutim)
			else :
				print(str(image_list_filter[0])+' is the only element. Pass.')
				pass
		elif large == True :
			if len(image_list_filter) > 1 : 
				for i in range(len(image_list_filter)):
					inim =  image_list_filter[i]
					part = inim.split('-')
					part[0] = 'Ref'
					newref = '-'.join(part)
					os.system('cp '+refim+' '+newref)
					images_to_align = newref
					ref_image = inim
					print('ref : '+inim) 
					outputshape	   = alipy.align.shape(ref_image)
					images_to_align_1 = []
					images		= images_to_align
					idfile = gregister([images])
					images_to_align_1.append(images)
					outim = idfile[0][:-5]+'_gregister.fits'
					newoutim = 'g'+idfile[0][:-5]+'.fits'
					os.system('cp ./alipy_out/'+outim+' ./'+newoutim)
			else :
				print(str(image_list_filter[0])+' is the only element. Pass.')
				pass
	os.system('mv alipy_out/*.fits .')
	os.system('rm -r alipy_out')
	print('Done. \a')
'''
#=================================================================
def wcsremap(srcim, tmpim, outim, verbose=False):
	'''
	1. Description
	run wcsremap-1.0.1 by Andrew Becker 
	wcsremap -template template.fits -source input.fits -outIm input_remapped.fits 

	srcim	   (Source)   : Larger  pixscale -> no change
	tmpim	   (Template) : Smaller pixscale -> change

	Small pixscale image ==> Large pixscale image
	'''
	import glob
	import os
	if verbose == True :
		com		= 'wcsremap -template '+tmpim+' -source '+srcim+' -outIm '+outim+' -v'
	elif verbose == False :
		com		= 'wcsremap -template '+tmpim+' -source '+srcim+' -outIm '+outim
	os.system(com)
	print(outim)
#-----------------------------------------------------------------
def run_wcsremap(imlist_name, refim, large=False):
	"""
	large = True : Sci pix is larger or the same with ref.
	large = False : Sci pix is smaller or when you make stacking image of one epoch.
	"""
	import glob
	import os
	imlist  = glob.glob(imlist_name)
	imlist.sort()
	if large == False :
		# pixscale of sciim is small than refim
		# sciim will be controlled.
		for i in range(len(imlist)):
			inim = imlist[i]
			parts	= inim.split('-')
			parts[0]= 'Remap'
			outim	= '-'.join(parts)	
			os.system('rm '+outim)
			wcsremap(inim, refim, outim, verbose=False)
	elif large == True :
		# pixscale of sciim is larger than refim
		# refim will be controlled.
		for i in range(len(imlist)):
			inim = imlist[i]
			parts	= inim.split('-')
			parts[0]= 'Remap'
			outim	= '-'.join(parts)
			if inim == outim :
				parts[0] = 'Ref_Remap'
				outim	= '-'.join(parts)
			os.system('rm '+outim)
			#os.system('cp '+refim+' '+outim)
			wcsremap(refim ,inim, outim, verbose=False)
	print('Done.')
#=================================================================
def run_wcsremap_epoch(imlist_name, large=False, sep=5.):
	"""
	large = True : Sci pix is larger or the same with ref.
	large = False : Sci pix is smaller or when you make stacking image of one epoch.
	"""
	import os, glob
	from lgpy.SameEpoch import SameEpoch
	res = SameEpoch(imlist_name, sep=sep)
	if large == False :
		# pixscale of sciim is small than refim
		# sciim will be controlled.
		for i in range(len(res)):
			imlist = list(res[i])
			for j in range(len(imlist)):
				inim = imlist[j]
				parts	= inim.split('-')
				parts[0]= 'Remap'
				outim	= '-'.join(parts)	
				os.system('rm '+outim)
				wcsremap(inim, imlist[0], outim, verbose=False)
	print('Done.')
	'''
	elif large == True :
		# pixscale of sciim is larger than refim
		# refim will be controlled.
		for i in range(len(imlist)):
			inim = imlist[i]
			parts	= inim.split('-')
			parts[0]= 'Remap'
			outim	= '-'.join(parts)
			if inim == outim :
				parts[0] = 'Ref_Remap'
				outim	= '-'.join(parts)
			os.system('rm '+outim)
			#os.system('cp '+refim+' '+outim)
			wcsremap(refim ,inim, outim, verbose=False)
	'''

#=================================================================
def imcopy(imlist_name, region) :
	"""
	"""
	import os, sys
	import glob
	from pyraf import iraf

	imlist = glob.glob(imlist_name)
	imlist.sort()
	for i in range(len(imlist)) :
		inim = imlist[i]
		outim = 't'+inim
		iraf.imcopy(inim+region, outim)
	print('Done.')
#=================================================================
def skysex(input_image, obs_ccd='SAO_STX16803', backsize=128, backfiltersize=5, mode='single'):
	"""
	1. Description 
	: Sky background subtraction code using SExtractor sky estimation algorithm. When using only one image, mode should be 'single'. When using more than two images, mode should be 'multi'.
	
	2. Usage
	>>> SAO_skysex() 

	3. History
	2019.12	 Created by G. Lim
	2019.01.25  Added to SAO_process.py by G. Lim
	"""
	import glob
	import os,sys
	import numpy as np
	from astropy.io import ascii
	from astropy.io import fits

	obscat    = ascii.read('/home/lim9/anaconda3/lib/python3.7/site-packages/lgpy/obs_spec.txt')
	obs_idx   = np.where(obs_ccd == obscat['obs_ccd'])[0]

	pixscale  = float(obscat['pixelscale'][obs_idx])
	gain      = float(obscat['gain'][obs_idx])
	fov       = float(obscat['fov_a'][obs_idx])

	#pixscale = 0.311 # SAO SBIG stx16803
	config = '/data1/code/sex.config/'
	DETECT_MINAREA  = '5'
	DETECT_THRESH   = '5'
	ANALYSIS_THRESH = '5'
	DEBLEND_NTHRESH = '32'
	DEBLEND_MINCONT = '0.005'
	BACK_SIZE	   = str(backsize)
	BACK_FILTERSIZE = str(backfiltersize)
	BACKPHOTO_TYPE  = 'GLOBAL'
	BACKPHOTO_THICK = '24' # default = 24

	print('DETECT_MINAREA	: '+DETECT_MINAREA)
	print('DETECT_THRESH	 : '+DETECT_THRESH)
	print('ANALYSIS_THRESH   : '+ANALYSIS_THRESH)
	print('DEBLEND_NTHRESH   : '+DEBLEND_NTHRESH)
	print('DEBLEND_MINCONT   : '+DEBLEND_MINCONT)
	print('BACK_SIZE		 : '+BACK_SIZE)
	print('BACK_FILTERSIZE   : '+BACK_FILTERSIZE) 
	print('BACKPHOTO_TYPE	: '+BACKPHOTO_TYPE)
	print('BACKPHOTO_THICK   : '+BACKPHOTO_THICK)

	### Single
	if mode == 'single' :
		#input_image = input('Input image? : ')
		inim = input_image
		print(inim+' is entered.')
		infilter = inim.split('-')[-2]
		data, hdr = fits.getdata(inim, header=True)
		#gain = hdr['EGAIN']
		check = 'BACKGROUND,-BACKGROUND'
		checkim = inim[:-5]+'-sky.fits,'+inim[:-5]+'-skysub.fits'
		sexcom='sex -c '+config+'skysub.sex '+inim+' -CATALOG_NAME '+inim[:-5]+'.skysub.cat -PARAMETERS_NAME '+config+'skysub.param -STARNNW_NAME '+config+'default.nnw -FILTER_NAME '+config+'default.conv -GAIN '+str(gain)+' -DETECT_THRESH '+DETECT_THRESH+' -ANALYSIS_THRESH '+ANALYSIS_THRESH+' -DETECT_MINAREA '+DETECT_MINAREA+ ' -PIXEL_SCALE ' +str(pixscale)+ ' -BACK_SIZE '+BACK_SIZE+' -BACK_FILTERSIZE '+BACK_FILTERSIZE+' -BACKPHOTO_TYPE '+BACKPHOTO_TYPE + ' -BACKPHOTO_THICK ' +BACKPHOTO_THICK+' -CHECKIMAGE_TYPE '+check+' -CHECKIMAGE_NAME '+checkim
		print(sexcom)
		os.system(sexcom)
		txt=open(inim[:-5]+'.skysub.txt','w+')
		txt.write('#image DETECT_MINAREA DETECT_THRESH ANALYSIS_THRESH DEBLEND_NTHRESH DEBLEND_MINCONT BACK_SIZE BACK_FILTERSIZE BACKPHOTO_TYPE BACKPHOTO_THICK'+'\n')
		txt.write(inim+' '+DETECT_MINAREA+' '+DETECT_THRESH+' '+ANALYSIS_THRESH+' '+DEBLEND_NTHRESH+' '+DEBLEND_MINCONT+' '+BACK_SIZE+' '+BACK_FILTERSIZE+' '+BACKPHOTO_TYPE+' '+BACKPHOTO_THICK+'\n')
		txt.close()
		print(checkim, ' are created.')

	### multi
	elif mode == 'multi' : 
		image = glob.glob('Cal*.fits')
		image.sort()
		os.system('rm ./Cal*sky*')
		for i in range(len(image)):
			inim = image[i]
			infilter = inim.split('-')[-2]
			data, hdr = fits.getdata(inim, header=True)
			#gain = hdr['EGAIN']
			#gain = 1.24
			check = 'BACKGROUND,-BACKGROUND'
			checkim = inim[:-5]+'-sky.fits,'+inim[:-5]+'-skysub.fits'

			sexcom='sex -c '+config+'skysub.sex '+inim+' -CATALOG_NAME '+inim[:-5]+'.skysub.cat -PARAMETERS_NAME '+config+'skysub.param -STARNNW_NAME '+config+'default.nnw -FILTER_NAME '+config+'default.conv -GAIN '+str(gain)+' -DETECT_THRESH '+DETECT_THRESH+' -ANALYSIS_THRESH '+ANALYSIS_THRESH+' -DETECT_MINAREA '+DETECT_MINAREA+ ' -PIXEL_SCALE ' +str(pixscale)+ ' -BACK_SIZE '+BACK_SIZE+' -BACK_FILTERSIZE '+BACK_FILTERSIZE+' -BACKPHOTO_TYPE '+BACKPHOTO_TYPE + ' -BACKPHOTO_THICK ' +BACKPHOTO_THICK+' -CHECKIMAGE_TYPE '+check+' -CHECKIMAGE_NAME '+checkim
			print(sexcom)
			os.system(sexcom)

			txt=open(inim[:-5]+'.skysub.txt','w+')
			txt.write('#image DETECT_MINAREA DETECT_THRESH ANALYSIS_THRESH DEBLEND_NTHRESH DEBLEND_MINCONT BACK_SIZE BACK_FILTERSIZE BACKPHOTO_TYPE BACKPHOTO_THICK'+'\n')
			txt.write(inim+' '+DETECT_MINAREA+' '+DETECT_THRESH+' '+ANALYSIS_THRESH+' '+DEBLEND_NTHRESH+' '+DEBLEND_MINCONT+' '+BACK_SIZE+' '+BACK_FILTERSIZE+' '+BACKPHOTO_TYPE+' '+BACKPHOTO_THICK+'\n')
			txt.close()
	print(checkim, ' are created.')
#=================================================================
def imscale(inim_list, zp0=30):
	"""
	Multiply constant to match zero point of images to combine taken in different sky, weather, exposure time.

	>>> imscale(imlist_name, refim=0)

	imlist_name : the image list to match zp.
	refim : the image to be reference. (0 = The first one, or specify filename.)
	"""
	import glob
	import os, sys
	import numpy as np
	from astropy.io import fits
	#coadd = glob.glob('coadd*.R.fits')
	coadd = glob.glob(inim_list)
	coadd.sort()
	ZP, ZP0 = [], zp0
	factor = []
	for i in range(len(coadd)) :
		inim = coadd[i]
		print(inim)
		data, hdr = fits.getdata(inim, header=True)
		ZP.append(hdr['ZP'])	
		factor_dum = 10.**(0.4*(ZP0 - ZP[i]))
		print(str(factor_dum))
		factor.append(factor_dum)
		newdata_name = 'f'+inim[:-5]+'.fits'
		newdata = data*factor_dum
		hdr['FLXSCALE'] = factor_dum 
		print('Writing new data as '+newdata_name+'...')
		fits.writeto(newdata_name, newdata, header=hdr, overwrite=True)
	print('Done.')
#=================================================================
def imcombine(imlist_name, wcsremap = False):
	"""
	1. Description
	: Stacking images in each object.
	(1) classify objects.
	(2) classify filters.
	(3) 
	2. Usage
	>>> SAO_imcombine()
	
	3. History
	2019.04.28 G. Lim added this code from separated code.
	"""
	import glob
	import os, sys
	import subprocess
	import numpy as np 
	from pyraf import iraf
	from astropy.io import fits

	lists= glob.glob(imlist_name)
	lists.sort()

	# (1) Object classification

	#lists = glob.glob('Cal*gre*.fits')
	objects = []
	for i in range(len(lists)) :
		inim = lists[i]
		split = inim.split('-')
		hdr = fits.getheader(inim)
		objects.append(hdr['object'])
	objectset = set(objects)
	objectlist = list(sorted(objectset))
	objectlist.sort()
	obj = 0
	for obj in range(len(objectlist)):
		object_name = objectlist[obj]
		if wcsremap == False :
			image_list = glob.glob('Cal*'+object_name+'*gre*.fits')
		elif wcsremap == True :
			image_list = glob.glob('Remap*'+object_name+'*0.fits')
		image_list.sort()
		k=0
		allfilter = []
		for k in range(len(image_list)) :
			hdr = fits.getheader(image_list[k])
			allfilter.append(hdr['filter'])
		filterset = set(allfilter)
		infilter = list(sorted(filterset))
		band = 0
		for band in range(len(infilter)):
			if wcsremap == False :
				list_name = 'Cal*'+object_name+'*'+infilter[band]+'*gre*.fits'
			elif wcsremap == True :
				list_name = 'Remap*'+object_name+'*'+infilter[band]+'*0.fits'
			image_list_gre = glob.glob(list_name)
			image_list_gre.sort()
			os.system('ls '+list_name+' > calibrated.list')
			hdr_gre = fits.getheader(image_list_gre[0])
			#mjd = []
			#for i in range(len(image_list_gre)):
			#	hdr_gre = fits.getheader(image_list_gre[i])
			#	mjd.append(hdr_gre['mjd'])
			#MJD = np.mean(mjd)	

			UTDATE = hdr_gre['date-obs'][0:10]
			UTSTART = hdr_gre['date-obs'][11:19]

			newimage = 'comCalib-SAO_STX16803-'+object_name+'-'+str(UTDATE[0:4])+str(UTDATE[5:7])+str(UTDATE[8:10])+'-'+str(UTSTART[0:2])+str(UTSTART[3:5])+str(UTSTART[6:8])+'-'+infilter[band]+'-'+str(int(hdr['exptime']*len(image_list_gre)))+'.fits'
			iraf.imcombine('@'+list_name, newimage, combine = 'median', reject='ccdclip', scale='none', zero='mode')

	print('Done.')
#=================================================================
def imcombine_set(imlist_name, reject='crreject'):
	"""
	Combine specific image set and make stacked image.
	> imcombine_set('DfCal*20190502-05*.fits')
	MJD will be added using median value of image set.

	none      - No rejection
	minmax    - Reject the nlow and nhigh pixels
	ccdclip   - Reject pixels using CCD noise parameters
	crreject  - Reject only positive pixels using CCD noise parameters
	sigclip   - Reject pixels using a sigma clipping algorithm
	avsigclip - Reject pixels using an averaged sigma clipping algorithm
	pclip     - Reject pixels using sigma based on percentiles
	"""
	import glob
	import os, sys
	import subprocess
	import numpy as np 
	from pyraf import iraf
	from astropy.io import fits
	from astropy.time import Time

	imlist = glob.glob(imlist_name)
	imlist.sort()
	mjd, exptime, date = [], [], []
	camera = imlist[0].split('-')[1]
	for i in range(len(imlist)) :
		inim = imlist[i]
		hdr  = fits.getheader(inim)
		mjd.append(hdr['MJD'])
		exptime.append(float(hdr['exptime']) )
		date.append(hdr['DATE-OBS'])
	
	combined_mjd	 = np.mean(mjd)	
	combined_exptime = np.sum(exptime)
	ncombine		 = len(imlist)
	part0 = imlist[0].split('-')[0]
	t= Time(combined_mjd, scale = 'utc', format='mjd')
	combined_date = t.fits[0:-4]

	print('MJD of combined image = '+str(combined_mjd))
	print('UT of combined image = '+str(combined_date))
	os.system('ls '+imlist_name+' > comb.list')
	object_name = hdr['object'] 
	if camera == 'LOAO' :
		band = hdr['filter'][0]
	else:
		band = hdr['filter']
	newhdr = hdr
	newhdr['MJD'] =  combined_mjd
	newhdr['EXPTIME'] = combined_exptime
	newhdr['NCOMBINE'] = ncombine
	try :
		newhdr['UTDATE']   = combined_date.split('T')[0]
		newhdr['UTSTART']  = combined_date.split('T')[1]
		newimage = 'com'+part0+'-'+camera+'-'+object_name+'-'+str(newhdr['UTDATE'][0:4])+str(newhdr['UTDATE'][5:7])+str(newhdr['UTDATE'][8:10])+'-'+str(newhdr['UTSTART'][0:2])+str(newhdr['UTSTART'][3:5])+str(newhdr['UTSTART'][6:8])+'-'+band+'-'+str(int(combined_exptime))+'.fits'
		print(newimage)
	except :
		newhdr['DATE-OBS'] = combined_date
		newimage = 'com'+part0+'-'+camera+'-'+object_name+'-'+str(newhdr['DATE-OBS'][0:4])+str(newhdr['DATE-OBS'][5:7])+str(newhdr['DATE-OBS'][8:10])+'-'+str(newhdr['DATE-OBS'][11:13])+str(newhdr['DATE-OBS'][14:16])+str(newhdr['DATE-OBS'][17:20])+'-'+band+'-'+str(int(combined_exptime))+'.fits'
		print(newimage)
	iraf.imcombine('@comb.list', newimage, combine = 'median', reject=reject, scale='none', zero='mode')
	newdata = fits.getdata(newimage)
	fits.writeto(newimage, newdata, header=newhdr, overwrite=True)
	print('New image and header is entered.')
#=================================================================
def imcombine_epoch(imlist_name, part0, sep=5.):
	"""
	Make images taken in the same epoch into one set and combine.
	You can control the separation. Default is 5 minutes.

	imlist_name = 'Remap*.fits'
	part0 = Calib, Remap,... front name of combined image
	sep = 5. minutes
	"""
	import glob
	import os, sys
	import numpy as np 
	from pyraf import iraf
	from lgpy.SameEpoch import SameEpoch
	from astropy.io import fits
	from astropy.time import Time
	from astropy.table import Table

	res = SameEpoch(imlist_name, sep=sep)
	for ii in range(len(res)) :
		image_to_com = list(res[ii])
		camera = image_to_com[0].split('-')[1]
		object_name = image_to_com[0].split('-')[2]
		band = image_to_com[0].split('-')[5]
		print(image_to_com)
		print('will be set as the same epoch.')
		mjd_to_com   = []
		exp_to_com   = []
		for iii in range(len(image_to_com)) :
			mjd_to_com.append(float(fits.getheader(image_to_com[iii])['mjd']))
			exp_to_com.append(float(fits.getheader(image_to_com[iii])['exptime']))
		new_exp     = np.sum(exp_to_com)
		mean_mjd	= np.mean(mjd_to_com)
		t			= Time(mean_mjd, format = 'mjd')
		new_dateobs = t.isot
		yy, mm, dd  = t.isot[0:4], t.isot[5:7], t.isot[8:10]
		hh, min, ss = t.isot[11:13], t.isot[14:16],  t.isot[17:19]

		newimage	 = 'com'+part0+'-'+camera+'-'+object_name+'-'+str(yy)+str(mm)+str(dd)+'-'+str(hh)+str(min)+str(ss)+'-'+band+'-'+str(int(new_exp))+'.fits'
		iraf.imcombine(','.join(image_to_com), newimage, combine = 'median', reject='avsigclip', scale='none', zero='mode')
#=================================================================
#------------------------------------------------------------
def getimages(ra,dec,size=240,filters="grizy"):
	"""Query ps1filenames.py service to get a list of images

	ra, dec = position in degrees
	size = image size in pixels (0.25 arcsec/pixel)
	filters = string with filters to include
	Returns a table with the results
	"""	
	from astropy.table import Table
	service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
	url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
		   "&filters={filters}").format(**locals())
	table = Table.read(url, format='ascii')
	return table
#------------------------------------------------------------
def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
	"""Get URL for images in the table


	ra, dec = position in degrees
	size = extracted image size in pixels (0.25 arcsec/pixel)
	output_size = output (display) image size in pixels (default = size).
				  output_size has no effect for fits format images.

	filters = string with filters to include
	format = data format (options are "jpg", "png" or "fits")
	color = if True, creates a color image (only for jpg or png format).
			Default is return a list of URLs for single-filter grayscale images.
	Returns a string with the URL
	"""
	import numpy as np
	import requests 
	from astropy.io.votable import parse_single_table 

	if color and format == "fits":
		raise ValueError("color images are available only for jpg or png formats")
	if format not in ("jpg","png","fits"):
		raise ValueError("format must be one of jpg, png, fits")

	table	= getimages(ra, dec, size=size, filters=filters)
	url		= (	"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
				"ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
	if output_size:
		url = url + "&output_size={}".format(output_size)
	# sort filters from red to blue
	flist = ["yzirg".find(x) for x in table['filter']]
	table = table[np.argsort(flist)]
	if color:
		if len(table) > 3:
			# pick 3 filters
			table = table[[0,len(table)//2,len(table)-1]]
		for i, param in enumerate(["red","green","blue"]):
			url = url + "&{}={}".format(param,table['filename'][i])
	else:
		urlbase = url + "&red="
		url = []
		for filename in table['filename']:
			url.append(urlbase+filename)
	return url
#------------------------------------------------------------
def ps1cut(imlist_name, infilter='r', fov=21.231):
	"""
	1. Description 
	: Sky background subtraction code using SExtractor sky estimation algorithm. When using only one image, mode should be 'single'. When using more than two images, mode should be 'multi'.
	
	2. Usage
	>>> SAO_skysex() 

	3. History
	2019.12	 Created by G. Lim
	2019.01.25  Added to SAO_process.py by G. Lim
	"""
	import glob
	import os, sys
	import numpy as np
	from astropy.io import fits
	from lgpy.hdrcheck import wcscenter
	from astropy.time import Time
	
	from pyraf import iraf
	#SAO_fov	= 21.231
	#MAO_fov = 18.3
	pixscale   = 0.25	# PS1
	ext_factor = 1.2

	imlist	 = glob.glob(imlist_name)
	imlist.sort()

	for i in range(len(imlist)):
		inim		   = imlist[i]		
		hdr			= fits.getheader(inim)
		obj			= hdr['object'] 
		#CRVAL1, CRVAL2 = hdr['CRVAL1'], hdr['CRVAL2']
		RA, DEC = wcscenter(inim)
		
		
		fov_pix		= ((fov*60.)/pixscale)*ext_factor
		if fov_pix > 6000 : 
			fov_pix = 6000
		param_geturl   = dict(  ra		  = RA,
								dec		 = DEC,
								size		= int(fov_pix),
								output_size = None,
								filters	 = infilter,
								format	  = "fits")
		#infilter	   = list(param_geturl.values())[3]
		try :
			url			= geturl(**param_geturl)
		except :
			try :
				url		= geturl(**param_geturl)
			except :
				try : 
					url		= geturl(**param_geturl)
				except :
					url		= geturl(**param_geturl)	
		save_dir   = './'+infilter
		os.system('mkdir '+save_dir)
		fh			 = fits.open(url[0])
		newname		= 'Ref-PS1-'+obj+'-'+infilter+'.fits'
		fh.writeto(save_dir+'/'+newname, overwrite=True)
		pan, panhd	 = fits.getdata(save_dir+'/'+newname, header=True)
		pan0		   = np.nan_to_num(pan)
		fits.writeto(save_dir+'/'+newname, pan0, panhd, overwrite=True)
		ps1im, ps1hdr = fits.getdata(save_dir+'/'+newname, header=True)
		iraf.hedit(save_dir+'/'+newname, fields='MJD', value = ps1hdr['MJD-OBS'], add='yes', verify='no')
		iraf.hedit(save_dir+'/'+newname, fields='OBJECT', value = obj, add='yes', verify='no')		
		
		t			= Time( ps1hdr['MJD-OBS'], format = 'mjd')
		new_dateobs = t.isot
		yy, mm, dd  = t.isot[0:4], t.isot[5:7], t.isot[8:10]
		hh, min, ss = t.isot[11:13], t.isot[14:16],  t.isot[17:19]
		iraf.hedit(save_dir+'/'+newname, fields='DATE-OBS', value =str(yy)+'-'+str(mm)+'-'+str(dd)+'T'+str(hh)+':'+str(min)+':'+str(ss), add='yes', verify='no')
	print('Please check the image at '+save_dir+'/'+newname)
#=================================================================
def gregister_inv(imlist_name, infilter='r'):
	"""
	Alipy gregister for subtraction using PS1
	from https://obswww.unige.ch/~tewes/alipy/index.html and its tutorial
	usage : run /data3/SAO1m/code/
	small pix --> large pix

	1. Description 
	: Alipy gregister for subtraction using archive data. (Ex. PS1). Reference images will be aligned to input_images. 

	from https://obswww.unige.ch/~tewes/alipy/index.html 
	
	2. Usage
	>>> SAO_gregister_inv('comCal*.fits', infilter='r') 

	3. History
	2019.12	 Created by G. Lim
	2019.01.25  Added to SAO_process.py by G. Lim
	"""
	import alipy
	import glob
	import os,sys
	import numpy as np
	from astropy.io import fits

	def gregister(ref_image, images):
		id = alipy.ident.run(ref_image, images, visu=False)
		print("%20s : %20s, flux ratio %.2f" % (id[0].ukn.name, id[0].trans, id[0].medfluxratio))
		alipy.align.irafalign(id[0].ukn.filepath, id[0].uknmatchstars, id[0].refmatchstars, shape=outputshape, makepng=False)

	imlist = glob.glob(imlist_name)
	imlist.sort()
	os.system('mkdir bad_align')
	for i in range(len(imlist)) :
		inim = imlist[i]
		hdr  = fits.getheader(inim)
		obj = hdr['object']
		ref = 'Ref-PS1-'+obj+'-'+infilter+'.fits'
		ref_path = '/data3/IMSNG/IMSNGgalaxies/refimg/SAO/'+infilter+'/'
		images_to_align = ref_path+ref #sys.argv[1] # Public data

		images = inim[:-5]+'.ref.fits'
		print('==================================')
		print('For '+inim+'...')
		print(ref+' is copied to '+images)
		print('==================================')
		os.system('cp '+images_to_align+' '+images)
		outputshape = alipy.align.shape(inim)
		try :
			gregister(inim, [images])
		except :
			os.system('mv '+inim+' '+images+' bad_align/')
			print(inim+' is not gregistered.')
			pass
	os.system('mv alipy_out/*.fits .')
	os.system('rm -r alipy_out')
	print('Done. \a')

#=================================================================
def hotpants_public(imlist_name) :
	"""
	1. Description 
	: This code aims to convolution and subtraction of already-gregistered reference images which are downloaded from other public archives. ex) PanStarrs DR1. These public images should be gregistered in advance by running gregister code using Alipy-based Python code.
	(1) Read Calibrated images.
	(2) Read reference images of each Calibrated images.
	(3) Running hotpants code (Becker's code) setting upper and lower limit of counts, which is investigated from original public image, and force convolution on template images with -c keyword. For upper & lower count limit of PS1 image, maximum count is more than 1.e+06 and minimum count is larger than -10000. 
	
	2. Usage
	>>> SAO_hotpants_public('Calib*.fits') :

	3. History
	2018.12.12 Created by G.Lim for Maidanak image and PS1 public data.
	2018.12.14 Test for SAO image and PS1 public data (on going)
	"""
	import glob, os
	# (1)
	objlist = glob.glob(imlist_name)
	objlist.sort()

	# (2)
	ref = []
	for i in range(len(objlist)):
		ref_name = objlist[i][:-5]+'.ref_gregister.fits'
		ref.append(ref_name)
	ref.sort()

	# (3)
	infile=objlist
	for n in range(len(infile)):
		outfile='hd'+infile[n]
		convfile='hc'+infile[n]
		com = 'hotpants -c t -n i -iu 2000000 -tu 2000000 -il -10000 -tl -10000 -v 0 -inim '+infile[n]+' -tmplim '+ref[n]+' -outim '+outfile+' -oci '+convfile
		print(infile[n])
		os.system(com)
	print('All done, check it out!')
#=================================================================
def hotpants(imlist_name, refim, scipsf=3, refpsf=2) :
	"""
	1. Description 
	: This code aims to convolution and subtraction of already-gregistered reference images which are downloaded from other public archives. ex) PanStarrs DR1. These public images should be gregistered in advance by running gregister code using Alipy-based Python code.
	(1) Read Calibrated images.
	(2) Read reference images of each Calibrated images.
	(3) Running hotpants code (Becker's code) setting upper and lower limit of counts, which is investigated from original public image, and force convolution on template images with -c keyword. For upper & lower count limit of PS1 image, maximum count is more than 1.e+06 and minimum count is larger than -10000. 
	
	imlist_name   : input science image list
	refim		 : Template image 
	sigmatch = True : In the case of Sigma_image > Sigma_template, for better subtraction, you may try this option. 

	2. Usage
	>>> SAO_hotpants('Calib*.fits') :

	3. History
	2018.12.12 Created by G.Lim for Maidanak image and PS1 public data.
	2018.12.14 Test for SAO image and PS1 public data (on going)
	"""
	import glob, os
	import numpy as np
	# (1)
	objlist = glob.glob(imlist_name)
	objlist.sort()

	# (2)
	ref = refim

	# (3)
	for n in range(len(objlist)):
		infile   = objlist[n]
		outfile  ='hd'+infile
		convfile ='hc'+infile
		print(infile[n])
		#FWHM = 2.355 sigma
		sigma_template = refpsf/ 2.355
		sigma_image	= scipsf / 2.355
		if sigma_image < sigma_template : # FWHM_sci is better 
			print('Science image will be convolved.')
			com	   = 'hotpants -c i -n i -iu 1000000 -tu 1000000 -il -1000 -tl -1000 -v 0 -inim '+infile+' -tmplim '+ref+' -outim '+outfile+' -oci '+convfile
		elif sigma_image > sigma_template : # FWHM_ref is better 
			print('Template image will be convolved.')
			com = 'hotpants -c t -n i -iu 1000000 -tu 1000000 -il -1000 -tl -1000 -v 0 -inim '+infile+' -tmplim '+ref+' -outim '+outfile+' -oci '+convfile
			#sigma_match	= np.sqrt(sigma_image**2 - sigma_template**2)
			#ngflag		 = ' -ng 3 6 '+ '%.3f'%(0.5*sigma_match) + ' 4 '+ '%.3f'%(sigma_match) +' 2 ' +'%.3f'%(2.0*sigma_match)
			#com			= 'hotpants -c t -n i -iu 70000 -tu 70000 -il -1000 -tl -1000 -v 1 -inim '+infile+' -tmplim '+ref+' -outim '+outfile+' -oci '+convfile+' '+ngflag
		print(com)
		os.system(com)
	print('All done, check it out!')
#=================================================================
def hotpants_other(sciim, refim, conv='t', scale='i', ngflag=False, refpsf='', scipsf='') :
	import glob, os
	import numpy as np
	from astropy.io import fits

	outfile  ='hd'+sciim
	convfile ='hc'+sciim
	#FWHM = 2.355 sigma
	if ngflag == False :
		com = 'hotpants -c '+conv+' -n '+scale+' -iu 65000 -tu 65000 -il -4000 -tl -4000 -v 0 -inim '+sciim+' -tmplim '+refim+' -outim '+outfile+' -oci '+convfile+' -ng 3 6 0.7 4 1.5 2 3.00'
	elif ngflag == True:
		sigma_template = refpsf/ 2.355
		sigma_image	= scipsf / 2.355
		#sigma_match = np.sqrt(sigma_image**2 - sigma_template**2)
		if sigma_image > sigma_template :
			sigma_match = np.sqrt(sigma_image**2 - sigma_template**2)
		elif sigma_image <= sigma_template :
			sigma_match = np.sqrt(sigma_template**2 - sigma_image**2)
		com = 'hotpants -c '+conv+' -n '+scale+' -iu 65000 -tu 65000 -il -4000 -tl -4000 -v 0 -inim '+sciim+' -tmplim '+refim+' -outim '+outfile+' -oci '+convfile+' -ng 3 6 '+str(0.5*sigma_match)+' 4 '+str(sigma_match)+' 2 '+str(2.*sigma_match)#+' -tg '+str(tgain)+' -tr '+str(trdnoise)+' -ig '+str(igain)+' -ir '+str(irdnoise)+' -fom v -r 10 -ft 20 -omi badmask.fits -ko 2 -ssig 3 -ks 2 -kfm 0.99 -allm -oki kernel.fits -okn -convvar' 
	print(com)
	os.system(com)

	print('All done, check it out!')	

'hotpants -c t -n i -iu 25000 -tu 25000 -il 0 -tl 0 -v 0 -inim comRemap-SAO_STX16803-NGC5350-20190625-12362-R-600.fits -tmplim Remap-SAO_STX16803-NGC5350-20190625-12362-R-600.fits -outim hdcomRemap-SAO_STX16803-NGC5350-20190625-12362-R-600.fits -oci hccomRemap-SAO_STX16803-NGC5350-20190625-12362-R-600.fits -ng 3 6 nan 4 nan 2 nan -tg 1.3600000143051147 -tr 9.0 -ig 1.3600000143051147 -ir 9.0 -fom v -r 10 -ft 20 -omi badmask.fits -ko 2 -ssig 3 -ks 2 -kfm 0.99 -allm -oki kernel.fits -okn -convvar'