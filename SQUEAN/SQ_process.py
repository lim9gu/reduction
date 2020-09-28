#run /data3/IMSNG/GW/GW190425/SQ_process.py
def SQ_biascom():
	"""
	1. Description 
	: 
	2. Usage 
	: 
	>>> SAO_biascom()

	3. History
	    Created by G.Lim.

	"""
	import glob
	import os, sys
	import numpy as np
	from pyraf import iraf
	from astropy.io import fits
	curdir = os.getcwd()
	curdate = curdir.split('/')[-1]
	iraf.noao()
	iraf.imred()
	iraf.ccdred()
	iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
	iraf.chdir('./bias')
	input_name = 'bias.list'
	output_name = curdate+'_zero.fits'
	#os.system('ls cal*bias.fit > '+input_name)
	calibrations = glob.glob('2019*.fits')
	f    = open(input_name,'w+')
	for i in xrange(len(calibrations)):
		hdr      = fits.getheader(calibrations[i])
		IMAGETYP = hdr['OBSTYPE']
		if IMAGETYP == 'BIAS' :	
			f.write(calibrations[i]+'\n')
	f.close()
	print('Zerocombine is running...')
	iraf.zerocombine(input='@'+input_name, output=output_name, combine='median', reject='minmax', process='no', scale='none', ccdtype='' )
	print('Output master '+output_name+' is created.')
	os.system('cp '+output_name+' ../')
	iraf.chdir('../')
	iraf.dir('.')
#=================================================================
def SQ_darkcom() :
	"""
	1. Description 
	:

	2. Usage
	: 

	3. History
	
	"""
	import glob
	import os, sys
	import numpy as np	
	from pyraf import iraf
	from astropy.io import fits
	curdir = os.getcwd()
	curdate = curdir.split('/')[-1]
	iraf.noao()
	iraf.imred()
	iraf.ccdred()
	iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
	iraf.chdir('./dark')
	dark = glob.glob('20*.fits')
	allexptime = []
	for i in xrange(len(dark)) :
		hdr = fits.getheader(dark[i])
		allexptime.append(hdr['exptime'])
	expset = set(allexptime)
	exptime = list(sorted(expset))
	i=0
	for i in xrange(len(exptime)) :
		print('Find images with exptime of '+str(exptime[i]))
		imlist = []
		for j in xrange(len(dark)) :
			hdr = fits.getheader(dark[j])
			if hdr['exptime'] == exptime[i] :
				imlist.append(dark[j])
			else :
				pass
		print(imlist)
		output_name = curdate+'_dark'+str(int(exptime[i]))+'.fits'
		input_name = output_name[:-5]+'.list'
		f=open(input_name,'w+')
		for k in xrange(len(imlist)) : 
			f.write(imlist[k]+'\n')
		f.close()
		print('Darkcombine is running...')
		iraf.imstat('@'+input_name)
		iraf.darkcombine(input='@'+input_name, output=output_name, combine='median', reject='minmax', process='no', scale='none', ccdtype='' )
		os.system('cp '+output_name+' ../')
	os.system('rm d*.list')
	iraf.chdir('../')
	iraf.dir('.')
	print('Output master '+output_name+' is created.')
#=================================================================
def SQ_flatcom(flattype='sky', process='bias', camera='SQUEAN') :
	"""
	1. Description 
	: 

	2. Usage
	: 
	3. History
	
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
	curdate = curdir.split('/')[-1]
	#flattype = sys.argv[1] # dome or sky
	#if flattype == 'dome_Hal' :
	#	iraf.chdir('./domeflat_Hal')
	#elif flattype == 'dome_LED' :
	#	iraf.chdir('./domeflat_LED')
	if flattype == 'dome' :
		iraf.chdir('./domeflat')
	elif flattype == 'sky' :
		iraf.chdir('./skyflat')
	flat = glob.glob('20*.fits')
	flat.sort()
	if process == 'bias' : # Bias subtraction : mainly skyflat
		input_name = 'flat.list'
		k=0
		f=open(input_name,'w+')
		for k in xrange(len(flat)) : 
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
		for i in xrange(len(flat)) :
			hdr = fits.getheader(flat[i])
			allexptime.append(hdr['exptime'])
		expset = set(allexptime)
		exptime = list(sorted(expset))
		i=0
		for i in xrange(len(exptime)) :
			print('Find images with exptime of '+str(exptime[i]))
			imlist = []
			j=0
			for j in xrange(len(flat)) :
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
			for k in xrange(len(flat)) : 
				f.write(flat[k]+'\n')
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
		calflat = glob.glob('z*.fits')
	elif process == 'dark' :
		calflat = glob.glob('d*.fits')
	allfilter = []
	i=0
	for i in xrange(len(calflat)) :
		hdr = fits.getheader(calflat[i])
		allfilter.append(hdr['filter'])
	filterset = set(allfilter)
	infilter = list(sorted(filterset))
	i=0
	for i in xrange(len(infilter)) :
		print('Find images with filter of '+str(infilter[i]))
		#if infilter[i] == 'Sii' :
		#	infilter[i] = 'SII'
		#elif infilter[i] == 'Oiii' :
		#	infilter[i] = 'OIII'
		#elif infilter[i] == 'ha' :
		#	infilter[i] = 'Ha'
		imlist = []
		for j in xrange(len(calflat)) :
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
		for k in xrange(len(imlist)) : 
			f.write(imlist[k]+'\n')
		f.close()
		#flatjoin = ",".join(imlist)
		output_name = input_name[:-5]+'.fits'
		iraf.flatcombine(input='@'+input_name, output=output_name, combine='median', reject='minmax', process='no', scale='mode', ccdtype='')
		#iraf.flatcombine(input=flatjoin, output=output_name, combine='median', reject='minmax', process='no', scale='mode', ccdtype='')
		print(output_name+' is created. Normalizing...')
		data, newhdr = fits.getdata(output_name, header=True)
		x = np.mean(data)
		nimage = data/x
		newflat_name = curdate+'_n'+str(infilter[i])+'flat.'+flattype+'.fits'
		fits.writeto(newflat_name, nimage, header=newhdr, overwrite=True)
		os.system('cp '+newflat_name+' ../')
	print('Normalised master flats are created.')
	iraf.imstat(images='*n?flat.'+flattype+'.fits')
	os.system('rm *.list ?flat.fits')
	iraf.chdir('../')
	iraf.dir('./')
#=================================================================
def SQ_objpre(flattype='sky') :
	"""
	1. Description
	: 
 
	2. Usage
	: 

	3. History

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
	obj_list = '20*_0*.fits'
	obj = glob.glob(obj_list)
	# dark subtraction
	allexptime = []
	i=0
	for i in xrange(len(obj)) :
		hdr = fits.getheader(obj[i])
		allexptime.append(hdr['exptime'])
	expset = set(allexptime)
	exptime = list(sorted(expset))
	i, j, k = 0, 0, 0
	for i in xrange(len(exptime)) :
		print('Find images with exptime of '+str(exptime[i]))
		imlist = []
		for j in xrange(len(obj)) :
			hdr = fits.getheader(obj[j])
			if hdr['exptime'] == exptime[i] :
				imlist.append(obj[j])
			else :
				pass
		print(imlist)
		imlist.sort()
		print('Creating object list for dark subtraction...')
		f = open("obj"+str(int(exptime[i]))+".list", 'w+')
		for im in xrange(len(imlist)) :
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
	dobj = glob.glob('d'+obj_list)
	dobj.sort()
	allfilter = []
	i=0
	for i in xrange(len(dobj)) :
		hdr = fits.getheader(dobj[i])
		allfilter.append(hdr['filter'])
	filterset = set(allfilter)
	infilter = list(sorted(filterset))
	i, j, k = 0, 0, 0
	for i in xrange(len(infilter)) :
		print('Find images with filter of '+str(infilter[i]))
		imlist = []
		imlist.sort()
		for j in xrange(len(dobj)) :
			hdr = fits.getheader(dobj[j])
			if hdr['filter'] == infilter[i] :
				imlist.append(dobj[j])
			else :
				pass
		print(imlist)

		g = open("obj"+infilter[i]+".list", 'w+')
		for im in xrange(len(imlist)) :
			g.write(imlist[im]+"\n")
		g.close()
		#dimjoin = ",".join(imlist)
		#doutimjoin = ",f".join(imlist)
		#doutimjoin = 'f'+doutimjoin
		print('Performing flat fielding...')
		nflats = glob.glob('*n'+str(infilter[i])+'flat.'+flattype+'.fits')[0]
		#nflats = 'n'+str(infilter[i])+'flat.'+flattype+'.fits'
		iraf.imarith(operand1='@obj'+infilter[i]+'.list', op='/', operand2=nflats, result='f@obj'+infilter[i]+'.list')
		#iraf.imarith(operand1=dimjoin, op='/', operand2=nflats, result=doutimjoin)
	print('Flat fielding is finished. Check the images.')
#=================================================================
def SQ_astrometry() :
	"""
	1. Description 
	: 

	2. Usage
	>>> 

	3. History

	"""
	import os,sys
	import glob
	import subprocess
	import numpy as np
	addlist = glob.glob('fd*.fits')
	addlist.sort()
	sexconfig = '/data3/SAO1m/code/astrom.config/astrometry.net.sex'
	print('Solving WCS using Astrometry.net...')
	for n in range(len(addlist)):
		com='solve-field '+addlist[n]+' --cpulimit 300 --overwrite --use-sextractor  --sextractor-config '+sexconfig+' --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_AUTO --sort-ascending --scale-unit arcsecperpix --scale-low 0.27 --scale-high 0.28 --no-remove-lines --uniformize 0 --no-plots  --new-fits a'+addlist[n]+' --temp-dir .\n'
		#com='solve-field '+addlist[n]+' --cpulimit 300 --overwrite --use-sextractor --scale-unit arcsecperpix --scale-low 0.27 --scale-high 0.28 --no-remove-lines --uniformize 0 --no-plots  --new-fits a'+addlist[n]+' --temp-dir .\n'
		print(com)
		print(str(n)+' th of '+str(len(addlist)))
		os.system(com) 
	orinum = subprocess.check_output('ls fd*.fits | wc -l', shell=True)
	resnum = subprocess.check_output('ls afd*.fits | wc -l', shell=True)
	print("from "+str(orinum[:-1])+" files , "+str(resnum[:-1])+" files are solved.")
	print("All done.")
	os.system('rm tmp*')
	os.system('rm *.wcs *.rdls *.corr *.xyls *.solved *.axy *.match ')
	print('Astrometry process is complete.')
#=================================================================
def SQ_scamp(imlist_name, refcat='USNO-B1') :
	"""
	# USNO-A1,USNO-A2,USNO-B1,
    # GSC-1.3,GSC-2.2,GSC-2.3,
    # TYCHO-2, UCAC-1,UCAC-2,UCAC-3,UCAC-4,
    # NOMAD-1, PPMX, CMC-14, 2MASS, DENIS-3,
    # SDSS-R3,SDSS-R5,SDSS-R6,SDSS-R7,
    # SDSS-R8, SDSS-R9
	"""
	import glob
	import os, sys
	from pyraf import iraf
	import subprocess
	import numpy as np
	import astropy.units as u
	from astropy.io import fits
	from astropy.io import ascii
	import astropy.coordinates as coord

	# SQUEAN
	#cd11   =    -7.66E-05 #/ Transformation matrix                          
	#cd12   =    -1.74E-06 #/ no comment                                     
	#cd21   =    -2.07E-06 #/ no comment                                     
	#cd22   =    7.68E-05 #/ no comment                                     

	cd11   = -7.9607532911290E-05 #/ Linear projection matrix                       
	cd12   =  -1.754777100142E-06 #/ Linear projection matrix                       
	cd21   =   -1.89497183878E-06 #/ Linear projection matrix                       
	cd22   =   7.969879099597E-05 #/ Linear projection matrix 

	sexconfig   = '/data3/SAO1m/code/astrom/astrom.sex'
	sexconv     = '/data3/SAO1m/code/astrom/astrom.conv'
	sexnnw      = '/data3/SAO1m/code/astrom/astrom.nnw'
	sexparam    = '/data3/SAO1m/code/astrom/astrom.param'
	psfconfig   = '/data3/SAO1m/code/psfex.config/prepsfex.sex'
	scampconfig = '/data3/SAO1m/code/astrom/astrom.scamp'

	#imlist      = glob.glob('Cal*gre*.fits')
	imlist      = glob.glob(imlist_name)
	imlist.sort()

	def wcsreset(inim,wcschar):
		iraf.wcsreset.image = inim
		iraf.wcsreset.wcs   = wcschar
		iraf.wcsreset.mode  = 'h'
		iraf.wcsreset()

	def hedit(inim,kword,value):
		iraf.hedit.images = inim
		iraf.hedit.fields = kword
		iraf.hedit.value  = value
		iraf.hedit.add    = 'yes'
		iraf.hedit.verify = 'no'
		iraf.hedit.update = 'yes'
		iraf.hedit.mode   = 'h'
		iraf.hedit()

	def scampastrom(inim):
		header = fits.getheader(inim)
		hdr    = header
		# ra    = header['ra']
		# dec   = header['dec']
		ra= '13:46:02.793'
		dec= '+55:42:03.29'
		rad   = coord.Angle(ra,unit=u.hour)	
		radd  = rad.degree
		decd  = coord.Angle(dec,unit=u.deg)
		decdd = decd.degree	
		# xpix  = hdr['NAXIS1']/2.
		# ypix  = hdr['NAXIS2']/2.
		xpix  = 685.65657
		ypix  = 347.35144
		threshold,minarea = 2,4
		output_cat = 'astromtest.cat'
		sexcom = 'sex -c '+sexconfig +' '+inim+' -PARAMETERS_NAME '+sexparam+' -DETECT_THRESH '+str(threshold) +' -DETECT_MINAREA '+str(minarea)+' -CATALOG_NAME '+output_cat+' -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+' -CATALOG_TYPE FITS_LDAC'
		scampcom='scamp -c '+scampconfig+' astromtest.cat '+ '-ASTREF_CATALOG '+refcat+' -CHECKPLOT_DEV PNG'

		# iraf setting
		wcsreset(inim,'physical')
		wcsreset(inim,'world')
		hedit(inim,'WAT0_001','system=image')
		hedit(inim,'WAT1_001', 'wtype=tan axtype=ra')
		hedit(inim,'WAT2_001', 'wtype=tan axtype=dec')
		hedit(inim,'RADECSYS', 'FK5')
		hedit(inim,'EQUINOX', '2000.')
		hedit(inim,'CTYPE1', 'RA---TAN')
		hedit(inim,'CTYPE2', 'DEC--TAN')
		hedit(inim,'CRVAL1', str(radd))
		hedit(inim,'CRVAL2', str(decdd))
		hedit(inim,'CRPIX1', str(xpix))
		hedit(inim,'CRPIX2', str(ypix))
		hedit(inim,'CD1_1', str(cd11))
		hedit(inim,'CD1_2', str(cd12))
		hedit(inim,'CD2_1', str(cd21))
		hedit(inim,'CD2_2', str(cd22))

		os.system(sexcom)
		os.system('ldactoasc '+output_cat+' > tmp.cat')
		tmpcat = ascii.read('tmp.cat')
		dolist = []
		nolist = []
		os.system('mkdir bad_align')
		if len(tmpcat['NUMBER']) < 5 : 
			print('Matched star number = ' + str(len(tmpcat['NUMBER'])))
			print('I will not touch this, it has stars less than 5. \n')
			nolist.append(inim+'\n')
			os.system('mv '+inim+' bad_ailgn/')
		else :
			print('Scamp will do this!! \n')
			dolist.append(inim+'\n')
			os.system(scampcom)
			hdr1 = hdr
			hdr1.extend(fits.Header.fromtextfile('astromtest.head'), update=True, update_first=True)
			fits.writeto('a'+inim,fits.getdata(inim),hdr1)
			print(inim,'astrometry work using SCAMP')

	for i in imlist : 
		os.system('delwcs '+i)
		os.system('rm a'+i)
		scampastrom(i)
#=================================================================
def SQ_gregister(imlist_name) :
	"""
	1. Description 
	: 
	
	2. Usage
	>>> SQ_gregister() 

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
	lists      = glob.glob(imlist_name)
	lists.sort()
	objects    = []
	for i in xrange(len(lists)) :
		hdr    = fits.getheader(lists[i])
		objects.append(hdr['object'])
	objectset  = set(objects)
	objectlist = list(sorted(objectset))
	
	def gregister(images):
		print('ref_image = ' + ref_image )
		print('Input_image ')
		id     = alipy.ident.run(ref_image, images, visu=False)
		print("%20s : %20s, flux ratio %.2f" % (id[0].ukn.name, id[0].trans, id[0].medfluxratio))
		alipy.align.irafalign(id[0].ukn.filepath, id[0].uknmatchstars, id[0].refmatchstars, shape=outputshape, makepng=False)

	ref_path = '/data3/IMSNG/GW/GW190425/SQUEAN/ps1/'

	### Filter classification
	obj = 0
	for obj in xrange(len(objectlist)):
		object_name = objectlist[obj]
		image_list  = glob.glob('Cal*'+object_name+'*0.fits')
		k           = 0
		allfilter   = []
		for k in xrange(len(image_list)) :
			hdr     = fits.getheader(image_list[k])
			allfilter.append(hdr['filter'])
		filterset   = set(allfilter)
		infilter    = list(sorted(filterset))	
		band        = 0
		for band in xrange(len(infilter)):
			image_list_filter     = glob.glob('Cal*'+object_name+'*'+infilter[band]+'*0.fits')
			if len(image_list_filter) > 1 : 
				images_to_align   = image_list_filter
				ref_image         = image_list_filter[1]
				outputshape       = alipy.align.shape(ref_image)
				images_to_align_1 = []
				n = 0
				for n in xrange(len(images_to_align)) :
					images        = images_to_align[n:n+1][0]
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
def SQ_fnamechange(camera='SQUEAN') :
	"""
	1. Description 
	: Change file name of WCS solving images (afd*.fits) using naming sequence following: Calib-SAO-CCD_name-OBJECT-UTDATE-UTSTART-FILTER-EXPTIME.fits. Then the code copies changed files to IMSNGgalaxies directory.

	2. Usage
	>>> SAO_fnamechange()

	3. History
	2018.03    Created by G.Lim.
	2018.12.21 Edited by G.Lim. Define SAO_fnamechange function.
	"""
	import glob
	import os, sys
	import numpy as np 
	import subprocess
	from astropy.io import fits
	lists = glob.glob('fd*.fits')
	lists.sort()
	for i in xrange(len(lists)):
		hdr = fits.getheader(lists[i])
		EXPTIME = int(hdr['exptime'])
		FILTER = hdr['filter']
		OBJECT = hdr['object']
		if (OBJECT == 'M51') | (OBJECT == 'M51a') :
			print('Object name is corrected.')
			OBJECT = 'M51A'
		UTDATE = hdr['DATE-OBS'][0:10]
		UTSTART = hdr['TIME-OBS'][0:8]
		newimage = 'Calib-MCD_'+camera+'-'+OBJECT+'-'+str(UTDATE[0:4])+str(UTDATE[5:7])+str(UTDATE[8:10])+'-'+str(UTSTART[0:2])+str(UTSTART[3:5])+str(UTSTART[6:8])+'-'+FILTER+'-'+str(EXPTIME)+'.fits'
		os.system('cp '+lists[i]+' '+newimage)
		print('Copy '+lists[i]+' to '+newimage+'.')
	print("Basic preprocessing of SQUEAN is finished.")
#=================================================================
def SQ_skysex(input_image, mode='single'):
	"""
	1. Description 
	: Sky background subtraction code using SExtractor sky estimation algorithm. When using only one image, mode should be 'single'. When using more than two images, mode should be 'multi'.
	
	2. Usage
	>>> SAO_skysex() 

	3. History
	2019.12     Created by G. Lim
	2019.01.25  Added to SAO_process.py by G. Lim
	"""
	import glob
	import os,sys
	import numpy as np
	from astropy.io import ascii
	from astropy.io import fits

	pixscale = 0.311 # SAO SBIG stx16803
	config = '/data3/SAO1m/code/sex.config/'

	DETECT_MINAREA  = '5'
	DETECT_THRESH   = '5'
	ANALYSIS_THRESH = '5'
	DEBLEND_NTHRESH = '32'
	DEBLEND_MINCONT = '0.005'
	BACK_SIZE       = '128'
	BACK_FILTERSIZE = '5'
	BACKPHOTO_TYPE  = 'GLOBAL'
	BACKPHOTO_THICK = '24' # default = 24

	print('DETECT_MINAREA    : '+DETECT_MINAREA)
	print('DETECT_THRESH     : '+DETECT_THRESH)
	print('ANALYSIS_THRESH   : '+ANALYSIS_THRESH)
	print('DEBLEND_NTHRESH   : '+DEBLEND_NTHRESH)
	print('DEBLEND_MINCONT   : '+DEBLEND_MINCONT)
	print('BACK_SIZE         : '+BACK_SIZE)
	print('BACK_FILTERSIZE   : '+BACK_FILTERSIZE) 
	print('BACKPHOTO_TYPE    : '+BACKPHOTO_TYPE)
	print('BACKPHOTO_THICK   : '+BACKPHOTO_THICK)

	### Single
	if mode == 'single' :
		#input_image = input('Input image? : ')
		inim = input_image
		print(inim+' is entered.')
		#infilter = inim.split('-')[-2]
		data, hdr = fits.getdata(inim, header=True)
		infilter = hdr['filter']
		gain = hdr['GAIN']
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
		image = glob.glob(input_image)
		os.system('rm ./Cal*sky*')
		for i in range(len(image)):
			inim = image[i]
			#infilter = inim.split('-')[-2]
			data, hdr = fits.getdata(inim, header=True)
			infilter = hdr['filter']
			gain = hdr['GAIN']
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
#========================================================
def SQ_imcombine():
	# Code for iraf imcombine (SAO 1m)
	### usage : run /data3/SAO1m/code/SAO_imcombine.py 

	import os, sys
	import numpy as np 
	from pyraf import iraf
	from astropy.io import fits
	import subprocess
	import glob
	'''
	objname = sys.argv[1]
	obsdate = sys.argv[2]
	infilter = sys.argv[3]

	infile='obj_'+objname+infilter+'.list'
	#subprocess.call('ls a*'+objname+'*'+infilter+'*.fit > '+infile,shell=True)
	subprocess.call('ls fd*'+objname+'*'+infilter+'*gregister.fits > '+infile,shell=True)
	#infile = 'obj.list'
	#subprocess.call('ls aCal*M101*R*gre*.fits > '+infile,shell=True)
	#outfile = 'Cal-M101-R-SNUO1m.imcombine.fits'
	files=np.genfromtxt(infile,usecols=(0),dtype=str)
	addlist=list(files)
	filenumber = len(addlist)
	outfile = 'Cal-'+obsdate+'-'+objname+'-'+infilter+'-SAO-'+str(filenumber)+'-imcombine.fits'
	print outfile
	iraf.imcombine('@'+infile, outfile, combine = 'median', reject='ccdclip', scale='none', zero='mode')
	#iraf.imcombine('@'+infile, outfile, combine = 'median', reject='minmax', scale='none', zero='mode', nlow=1., nhigh=3., lsigma=2., hsigma=2.)

	'''
	### Object classification
	#lists = glob.glob('fd*gre*.fits')
	lists = glob.glob('Cal*gre*.fits')
	objects = []
	objects.sort()
	i=0
	for i in xrange(len(lists)) :
		hdr = fits.getheader(lists[i])
		objects.append(hdr['object'])
	objectset = set(objects)
	objectlist = list(sorted(objectset))
	objectlist.sort()
	obj = 0
	for obj in xrange(len(objectlist)):
		object_name = objectlist[obj]
		#image_list = glob.glob('fd'+object_name+'*gre*.fits')
		image_list = glob.glob('Cal*'+object_name+'*gre*.fits')
		image_list.sort()
		k=0
		allfilter = []
		for k in xrange(len(image_list)) :
			hdr = fits.getheader(image_list[k])
			allfilter.append(hdr['filter'])
		filterset = set(allfilter)
		infilter = list(sorted(filterset))
		band = 0
		for band in xrange(len(infilter)):
			list_name = 'Cal*'+object_name+'*'+infilter[band]+'*gre*.fits'
			image_list_gre = glob.glob(list_name)
			image_list_gre.sort()
			os.system('ls '+list_name+' > calibrated.list')
			hdr_gre = fits.getheader(image_list_gre[0])
			UTDATE = hdr_gre['DATE-OBS'][0:10]
			UTSTART = hdr_gre['TIME-OBS'][0:8]
		
			newimage = 'Calib-McD_SQUEAN-'+object_name+'-'+str(UTDATE[0:4])+str(UTDATE[5:7])+str(UTDATE[8:10])+'-'+str(UTSTART[0:2])+str(UTSTART[3:5])+str(UTSTART[6:8])+'-'+infilter[band]+'-'+str(int(hdr['exptime']*len(image_list_gre)))+'.imcomb.fits'
			iraf.imcombine('@'+list_name, newimage, combine = 'median', reject='ccdclip', scale='none', zero='mode')

	print('Done.')
#=================================================================
def SQ_gregister_inv(imlist_name, infilter='i'):
	"""
	Alipy gregister for subtraction using PS1
	from https://obswww.unige.ch/~tewes/alipy/index.html and its tutorial
	usage : run /data3/SAO1m/code/
	small pix --> large pix

	1. Description 
	: Alipy gregister for subtraction using archive data. (Ex. PS1). Reference images will be aligned to input_images. 

	from https://obswww.unige.ch/~tewes/alipy/index.html 
	
	2. Usage
	>>> SAO_gregister_inv('Cal*imcomb.fits', infilter='r') 

	3. History
	2019.12     Created by G. Lim
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
	for i in xrange(len(imlist)) :
		inim = imlist[i]
		hdr  = fits.getheader(inim)
		obj = hdr['object']
		#ref = 'Ref-PS1-'+obj+'-'+infilter+'.fits'
		#ref_path = '/data3/IMSNG/IMSNGgalaxies/refimg/SAO/'+infilter+'/'
		#images_to_align = ref_path+ref #sys.argv[1] # Public data
		images_to_align = 'frame-i-003900-2-0759.fits'
		images = inim[:-5]+'.ref.fits'
		#print('==================================')
		#print('For '+inim+'...')
		#print(ref+' is copied to '+images)
		#print('==================================')
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
def SQ_hotpants_public(imlist_name) :
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
	#ref = []
	#for i in xrange(len(objlist)):
	#	ref_name = objlist[i][:-5]+'.ref_gregister.fits'
	#	ref.append(ref_name)
	#ref.sort()
	ref = ['hdRemap-SQUEAN-PGC58513-20190427-055032-i-180.fits']
	# (3)
	infile=objlist
	for n in xrange(len(infile)):
		outfile='hd'+infile[n]
		convfile='hc'+infile[n]
		com = 'hotpants -c t -n t -iu 2000000 -tu 2000000 -il -10000 -tl -10000 -v 0 -inim '+infile[n]+' -tmplim '+ref[n]+' -outim '+outfile+' -oci '+convfile
		print(infile[n])
		os.system(com)
	print('All done, check it out!')
#=================================================================
def retrieve_SQ(input_images, catalog='sdss', pixscale=0.276):
	"""
	1. Description 
	: Retrieve MAO reference catalog for MAO image. 
	(1) Do photometry roughly for extracting the RA, DEC in the image center.
	(2) 
		(i) Query SDSS dr12 stars with SQL command for clean photometry (stars).
		http://skyserver.sdss.org/dr12/en/help/docs/realquery.aspx#cleanStars
	    Clean stars have 17 < r_SDSS < 22, flagged with flag = 0, psfmagerr < 0.1, making them table catalog. For g, r, and i are close to AB. u_AB = u_SDSS - 0.04, z_AB = z_SDSS + 0.02.
		(ii) Query APASS9 stars using vizier system.

	(3) Convert SDSS mag [AB] into Johnson mag [Vega] using Lupton+05 formula below.
	For stars
	* used equation.(gri based, low sigma included equation is recommanded first.)
	   B = u - 0.8116*(u - g) + 0.1313;  sigma = 0.0095
	   *B = g + 0.3130*(g - r) + 0.2271;  sigma = 0.0107
	Bmag = gmag + 0.3130*(gmag - rmag) + 0.2271	# sigma = 0.0107
	   V = g - 0.2906*(u - g) + 0.0885;  sigma = 0.0129
	   *V = g - 0.5784*(g - r) - 0.0038;  sigma = 0.0054
	Vmag = gmag - 0.5784*(gmag - rmag) - 0.0038 # sigma = 0.0054
	   R = r - 0.1837*(g - r) - 0.0971;  sigma = 0.0106
	   *R = r - 0.2936*(r - i) - 0.1439;  sigma = 0.0072
	Rmag = rmag - 0.2936*(rmag - imag) - 0.1439 # sigma = 0.0072
	   *I = r - 1.2444*(r - i) - 0.3820;  sigma = 0.0078
	   I = i - 0.3780*(i - z)  -0.3974;  sigma = 0.0063
	Imag = rmag - 1.2444*(rmag - imag) - 0.3820 # sigma = 0.0078

	(i) Vega --> AB conversion (Frei & Gunn 1994, AJ, 108. 1476)
	#B = B(AB) + 0.163 +/- 0.004
	#V = V(AB) + 0.044 +/- 0.004
	#R = R(AB) - 0.055 +/- INDEF 
	 Conversion from AB magnitudes to Vega magnitudes:
    The following formulae convert between the AB magnitude systems and those based on Alpha Lyra:
        * V	=   V(AB) + 0.044	(+/- 0.004)
        * B	=   B(AB) + 0.163	(+/- 0.004)
        Bj	=  Bj(AB) + 0.139	(+/- INDEF)
        * R	=   R(AB) - 0.055	(+/- INDEF)
        * I	=   I(AB) - 0.309	(+/- INDEF)
         g	=   g(AB) + 0.013	(+/- 0.002)
         r	=   r(AB) + 0.226	(+/- 0.003)
         i	=   i(AB) + 0.296	(+/- 0.005)
         u'	=  u'(AB) + 0.0	        
         g'	=  g'(AB) + 0.0	        
         r'	=  r'(AB) + 0.0	        
         i'	=  i'(AB) + 0.0	        
         z'	=  z'(AB) + 0.0	        
        Rc	=  Rc(AB) - 0.117	(+/- 0.006)
        Ic	=  Ic(AB) - 0.342	(+/- 0.008)
    Source: Frei & Gunn 1994, AJ, 108, 1476 (their Table 2).

	(ii) Vega - AB Magnitude Conversion (Blanton+07)
		 U : m_AB - m_Vega =  0.79
		 B : m_AB - m_Vega = -0.09
		 V : m_AB - m_Vega =  0.02
		 R : m_AB - m_Vega =  0.21
		 I : m_AB - m_Vega =  0.45
		 J : m_AB - m_Vega =  0.91
		 H : m_AB - m_Vega =  1.39
		 K : m_AB - m_Vega =  1.85

	Used (ii) Blanton+07 equation since this is recent one.

	(4) Making table and save this into .cat file.

	2. Usage
	>>> run /data2/dwarf/code/candidate_selection/retrieve_mao.py 
	>>> retrieve_MAO(inim, pixscale)
		ex. retrieve_MAO('Calib*.fits', '0.267') 
		After running the code, you will see the saying that "Check the output file. Add '#' symbol in the first line." Follow this sign.

	3. History
	2019.03.13 Created by G.Lim based on retrive_sdss.py 
	"""
	import glob
	import os,sys
	import numpy as np
	from astropy.io import ascii
	from astroquery.sdss import SDSS
	from astropy.table import Table, Column

	import astropy.units as u
	from astropy.table import Table
	import astropy.coordinates as coord
	from astroquery.vizier import Vizier 
	from astropy.coordinates import Angle

	# (1) 
	inimlist   = glob.glob(input_images)
	apersize   = '3,5,7,9,11,13,15,17,19,21'
	detecthred = '5'
	analthred  = '5'
	config     = '/data3/SAO1m/code/sex.config/'
	for i in xrange(len(inimlist)):
		inim        = inimlist[i] 
		print(inim)
		precat_name = inim[:-5]+'.pre.cat'
		sexcom      = 'sex -c '+config+'default.sex '+inim+' -CATALOG_NAME '+precat_name+' -PARAMETERS_NAME '+config+'default.param -STARNNW_NAME '+config+'default.nnw -FILTER_NAME '+config+'default.conv -seeing_fwhm '+str(1.0)+' -DETECT_THRESH '+detecthred+' -ANALYSIS_THRESH '+analthred+' -PHOT_APERTURES '+apersize+' -PIXEL_SCALE '+str(pixscale)
		#sexcom      = 'sex -c '+config+'default.sex '+inim+' -CATALOG_NAME '+precat_name+' -PARAMETERS_NAME '+config+'default.param -STARNNW_NAME '+config+'default.nnw -FILTER_NAME '+config+'default.conv -seeing_fwhm '+str(1.0)+' -DETECT_THRESH '+detecthred+' -ANALYSIS_THRESH '+analthred+' -PIXEL_SCALE '+str(pixscale)
		print(sexcom)
		os.system(sexcom)

    # (2) SDSS dr12 or APASS
		if catalog == 'sdss' :
			incat           = ascii.read(precat_name)
			RA_arr, DEC_arr = incat['ALPHA_J2000'], incat['DELTA_J2000']
			RA, DEC         = np.median(RA_arr), np.median(DEC_arr)
			print('Central RA, DEC = ('+str(RA)+', '+str(DEC)+')')

			query = "SELECT top 10000 u,g,r,i,z,Err_u,Err_g,Err_r,Err_i, Err_z, ra,dec,flags_r from star where ra between "+str(np.min(RA_arr))+" and "+str(np.max(RA_arr))+" and dec between "+str(np.min(DEC_arr))+" and "+str(np.max(DEC_arr))+" and r > 17 and r < 22 and ((flags_r & 0x10000000)!=0) and ((flags_r & 0x8100000c00a4) = 0) and (((flags_r & 0x400000000000) = 0) or (psfmagerr_r <= 0.1)) AND (((flags_r & 0x100000000000) = 0) or (flags_r & 0x1000) = 0)"	
			print(query)
			sdss = SDSS.query_sql(query, data_release = 12)
			print(sdss[:5])
			### r_sdss -> Johnson_R Conversion (Robert Lupton 2007)
			sdss.colnames
			gmag, rmag, imag = sdss['g'],     sdss['r'],     sdss['i']
			gerr, rerr, ierr = sdss['Err_g'], sdss['Err_r'], sdss['Err_i']

			Bmag      = gmag + 0.3130*(gmag - rmag) + 0.2271	# sigma = 0.0107
			Berr      = np.sqrt(gerr**2 + 0.3130*(rerr**2+gerr**2) + 0.0107**2) 
			Bmag_AB   = Bmag - 0.09

			Vmag      = gmag - 0.5784*(gmag - rmag) - 0.0038 # sigma = 0.0054
			Verr      = np.sqrt(gerr**2 + 0.5784*(gerr**2+rerr**2) + 0.0054**2)
			Vmag_AB   = Vmag + 0.02

			Rmag      = rmag - 0.2936*(rmag - imag) - 0.1439 # sigma = 0.0072
			Rerr      = np.sqrt(rerr**2 + 0.2936*(rerr**2+ierr**2) + 0.0072**2)
			Rmag_AB   = Rmag + 0.21

			Imag = rmag - 1.2444*(rmag - imag) - 0.3820 # sigma = 0.0078
			Ierr = np.sqrt(rerr**2 + 1.2444*(rerr**2+ierr**2) + 0.0078**2)
			Imag_AB   = Imag + 0.45

			Bm = Column(name='Bmag', data=Bmag_AB)
			Be = Column(name='Berr', data=Berr)
			Vm = Column(name='Vmag', data=Vmag_AB)
			Ve = Column(name='Verr', data=Verr)
			Rm = Column(name='Rmag', data=Rmag_AB)
			Re = Column(name='Rerr', data=Rerr)
			Im = Column(name='Imag', data=Imag_AB)
			Ie = Column(name='Ierr', data=Ierr)

			sdss.add_column(Bm)
			sdss.add_columns([Be,Vm])
			sdss.add_columns([Ve,Rm,Re,Im,Ie])

			sdss.write(inim[:-5]+'.sdss.conv.cat', format='ascii')
			print("Check the output file. Add '#' symbol in the first line.")
		elif catalog == 'apass' :

			#input_precat = sys.argv[1]
			radius      = 1. # deg
			name        = precat_name.split('.')[2]
			precat      = ascii.read(precat_name)
			radeg       = precat['ALPHA_J2000']
			dedeg       = precat['DELTA_J2000']
			racen       = np.median(radeg)
			decen       = np.median(dedeg)
			FWHM_IMAGE  = precat['FWHM_IMAGE']

			comment = 'NAME'+'\t'+': '+name+'\n' \
					+ 'RA'+'\t'+': '+str(round(racen, 3))+'\n' \
					+ 'Dec'+'\t'+': '+str(round(decen, 3))+'\n' \
					+ 'Radius'+'\t'+': '+str(radius)+' deg'+'\n'*2 \
					+ 'LOADING APASS Catalog ...'+'\n'
			print(comment)
			outname = inim[:-5]+'.apass.conv.cat'
			Vizier.ROW_LIMIT    = -1
			query    = Vizier.query_region(coord.SkyCoord(ra=racen, dec=decen, \
										unit=(u.deg, u.deg), frame='icrs'), \
										width=str(radius*60)+'m', catalog=["APASS9"])
			dum      = query[0]
			colnames = dum.colnames
			for col in colnames:
				indx    = np.where( dum[col].mask == False )
				dum     = dum[indx]
			#   Vega    : B, V
			#   AB      : g, r, i
			#   Vega - AB Magnitude Conversion (Blanton+07)
			#   U       : m_AB - m_Vega = 0.79
			#   B       : m_AB - m_Vega =-0.09
			#   V       : m_AB - m_Vega = 0.02
			#   R       : m_AB - m_Vega = 0.21
			#   I       : m_AB - m_Vega = 0.45
			#   J       : m_AB - m_Vega = 0.91
			#   H       : m_AB - m_Vega = 1.39
			#   K       : m_AB - m_Vega = 1.85
			apasscat			= Table()
			apasscat['NUMBER']  = dum['recno']
			apasscat['RA_ICRS'] = dum['RAJ2000']
			apasscat['DE_ICRS'] = dum['DEJ2000']
			apasscat['Numb_obs']= dum['nobs']
			apasscat['Numb_img']= dum['mobs']
			apasscat['B-V']     = dum['B-V']    + (-0.09 - 0.02)
			apasscat['e_B-V']   = dum['e_B-V']
			apasscat['Bmag']    = dum['Bmag']   - 0.09  # [Vega] to [AB]
			apasscat['e_Bmag']  = dum['e_Bmag']
			apasscat['Vmag']    = dum['Vmag']   + 0.02  # [Vega] to [AB]
			apasscat['e_Vmag']  = dum['e_Vmag']
			apasscat['g']    = dum['g_mag']
			apasscat['e_g']  = dum['e_g_mag']
			apasscat['r']    = dum['r_mag']
			apasscat['e_r']  = dum['e_r_mag']
			apasscat['i']    = dum['i_mag']
			apasscat['e_i']  = dum['e_i_mag']

			Rmag_Vega           = apasscat['r'] - 0.2936*(apasscat['r'] - apasscat['i']) - 0.1439 # sigma = 0.0072
			Rerr      = np.sqrt(apasscat['e_r']**2 + 0.2936*(apasscat['e_r']**2+apasscat['e_i']**2) + 0.0072**2)
			Rmag_AB             = Rmag_Vega + 0.21 # [Vega] to [AB]

			Imag_Vega           = apasscat['r'] - 1.2444*(apasscat['r'] - apasscat['i']) - 0.3820 # sigma = 0.0078
			Ierr = np.sqrt(apasscat['e_r']**2 + 1.2444*(apasscat['e_r']**2+apasscat['e_i']**2) + 0.0078**2)
			Imag_AB             = Imag_Vega + 0.45 # [Vega] to [AB]
			apasscat['Rmag']    = Rmag_AB
			apasscat['e_Rmag']  = Rerr
			apasscat['Imag']    = Imag_AB
			apasscat['e_Imag']  = Ierr
			apasscat.rename_column('NUMBER', '#NUMBER_SE')
			apasscat.write(outname, format='ascii', overwrite=True)
			#print "Check the output file. Add '#' symbol in the first line." 

'''
import glob
imlist = glob.glob('Cal*.fits')
imlist.sort()
for i in xrange(len(imlist)):
	print imlist[i]
	retrieve_SQ(imlist[i], str(0.276))
print 'Done.'
'''

#=================================================================
#def SQ_PS1

#=================================================================
def SQ_fwhmzp_sdss12(inim, lowermag = 17, uppermag = 22, pixscale = 0.267, hist=False) :
	"""
	1. Description 

	2. Usage

	3. History
	"""
	import glob
	import os,sys
	import numpy as np
	from numpy import median
	from astropy.io import fits
	from astropy.io import ascii
	from termcolor import colored
	import matplotlib.pyplot as plt
	from astropy.stats import sigma_clip

	# (1)
	config     = '/data2/dwarf/code/sex.config/'
	apersize   = '3,5,7,9,11,13,15,17,19,21'
	detecthred = '5'
	analthred  = '5'

	#inim = sys.argv[1]
	if 'fscale' not in inim.split('.') :
		wgtim      = './'+inim[:-5]+'.weight.fits'
		refcat     = inim[:-5]+'.sdss.conv.cat'
	elif 'fscale' in inim.split('.') :
		wgtim      = './'+inim[:-12]+'.weight.fits'
		refcat     = inim[:-12]+'.sdss.conv.cat'
	print('Find weight images :'+wgtim)
	band           = inim.split('.')[3]
	infilter       = band+'mag'
	print('For '+infilter+'.')
	galaxy         = inim.split('.')[2]
	
	data, hdr      = fits.getdata(inim, header=True)
	exptime        = hdr['exptime']
	NAXIS1, NAXIS2 = hdr['NAXIS1'], hdr['NAXIS2']

	X_cen , Y_cen, XY_cut = NAXIS1/2, NAXIS2/2, 0.3
	X_off1, X_off2        = X_cen + NAXIS1*XY_cut, X_cen - NAXIS1*XY_cut 
	Y_off1, Y_off2        = Y_cen + NAXIS2*XY_cut, Y_cen - NAXIS2*XY_cut  

	precat_name = inim[:-5]+'.dum.cat'
	precom = 'sex -c '+config+'default.sex '+inim+' -CATALOG_NAME '+precat_name+' -PARAMETERS_NAME '+config+'default.param -STARNNW_NAME '+config+'default.nnw -FILTER_NAME '+config+'default.conv -GAIN '+str(gain)+' -seeing_fwhm '+str(1.0)+' -DETECT_THRESH '+detecthred+' -ANALYSIS_THRESH '+analthred+' -PHOT_APERTURES '+apersize+' -PIXEL_SCALE '+str(pixscale)+' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+wgtim+' -BACK_TYPE AUTO -BACKPHOTO_TYPE LOCAL -BACKPHOTO_THICK 24'
	print(precom)
	os.system(precom)
	precat        = ascii.read(precat_name)
	preX, preY    = precat['X_IMAGE'], precat['Y_IMAGE']
	premagauto    = precat['MAG_AUTO']
	premagautoerr = precat['MAGERR_AUTO']
	prefwhm       = precat['FWHM_IMAGE']*float(pixscale)
	preflag       = precat['FLAGS']
	prestell      = precat['CLASS_STAR']
	preindex = np.where((preX <= X_off1) & (preX >= X_off2) & (preY <= Y_off1) & (preY >= Y_off2) & (premagauto > -16) & (premagauto < -10) & (preflag == np.min(preflag)) )[0]
	preseeing_dum = np.median(prefwhm[preindex])

	#(2)
	zpcat    = inim[:-5]+".zp.cat"
	sexcom   = 'sex -c '+config+'default.sex '+inim+' -CATALOG_NAME '+zpcat+' -PARAMETERS_NAME '+config+'default.param -STARNNW_NAME '+config+'default.nnw -FILTER_NAME '+config+'default.conv -GAIN '+str(gain)+' -seeing_fwhm '+str(preseeing_dum)+' -DETECT_THRESH '+detecthred+' -ANALYSIS_THRESH '+analthred+' -PHOT_APERTURES '+apersize+' -PIXEL_SCALE '+str(pixscale)+' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+wgtim+' -BACK_TYPE AUTO -BACKPHOTO_TYPE LOCAL -BACKPHOTO_THICK 24'
	print(sexcom)
	os.system(sexcom)
	mergecat = inim[:-5]+".merge.sdss.cat"
	matchcom = "/home/lim9/stilts/stilts tskymatch2 ifmt1=ascii ifmt2=ascii in1="+refcat+" in2="+zpcat+" out="+mergecat+" ra1=ra dec1=dec ra2=col12 dec2=col13 error=2 join=1and2 ofmt=ascii omode=out" # write '#' in the first line in refcat.
	print(matchcom)
	os.system(matchcom)
	try : 
		incat = ascii.read(mergecat)
	except :
		print(colored('merge.cat is not created. Check your column number of refcat or photcat!', 'red'))
		pass
	X          = incat['col10']          # Cut X image coordinate to include image boundary. 
	Y          = incat['col11']          # Cut Y image coordinate to include image boundary.
	nomadmag   = incat[infilter]         # Johnson [AB] converted from SDSS.	
	nomadmage  = incat[band+'err']       # Johnson [AB] mag error propagated
	magauto    = incat['col4']           # Instrumental mag in input image.
	magautoerr = incat['col5']           # Instru. mag error in input image.
	fwhm       = incat['col16']*float(pixscale) # FWHM image gaussian core [arcsec]
	flag       = incat['col15']          # Flag = 0 : No problem
	stell      = incat['col17']      

	# Star selection criteria.
	try :
		index1 = np.where( (X >= X_off2) & (X <= X_off1) & (Y >= Y_off2) & (Y <= Y_off1) )[0]
		index2 = np.where( (X >= X_off2) & (X <= X_off1) & (Y >= Y_off2) & (Y <= Y_off1) & (flag == np.min(flag)) )[0]
		index3 = np.where( (X >= X_off2) & (X <= X_off1) & (Y >= Y_off2) & (Y <= Y_off1) & (flag == np.min(flag)) & (magautoerr < 0.1) & (nomadmage < 0.05))[0]
		index4 = np.where( (X >= X_off2) & (X <= X_off1) & (Y >= Y_off2) & (Y <= Y_off1) & (flag == np.min(flag))  & (magautoerr < 0.1) & (nomadmage < 0.05) & (nomadmag > lowermag) & (nomadmag < uppermag))[0]
		print(str(len(incat[index1]))+' stars are selected.')
		print(str(len(incat[index2]))+' stars are selected.')
		print(str(len(incat[index3]))+' stars are selected.')
		print(str(len(incat[index4]))+' stars are finally selected for measuring FWHM & ZP.')
	except :
		print(colored(inim + ' No stars within criteria.', 'red'))
		pass

	# Zero point value.
	zplist          = nomadmag[index4] - magauto[index4]
	sigma           = 2
	zplist_clip     = sigma_clip(zplist , sigma=sigma, iters=None, cenfunc=median, copy=False)
	indx_alive      = np.where( zplist_clip.mask == False )[0]
	indx_exile      = np.where( zplist_clip.mask == True )[0]
	intbl_alive     = incat[index4][indx_alive]
	intbl_exile     = incat[index4][indx_exile]
	print(str(len(intbl_alive))+' stars are within '+str(sigma)+' sigma.')
	print(str(len(intbl_exile))+' stars are clipped.')
	zplist_alive    = intbl_alive[infilter] - intbl_alive['col4']
	zp              = np.median(zplist_alive)
	zper            = np.std(zplist_alive)
	seeing_dum      = np.median(intbl_alive['col16']*float(pixscale))

	print('ZP median        = '+str(round(zp,3)))
	print('ZP stddev        = '+str(round(zper,3)) )
	print('Used star number = '+str(len(incat[indx_alive])))
	print('Seeing ["]       = '+str(round(seeing_dum,3)))

	hdr['Seeing']  = (round(seeing_dum,3),'FWHM in arcsec estimated from SExtractor')
	hdr['ZP']      = (round(zp,3),'Photometric zero-point [mag] using APASS catalog')
	hdr['ZPstd']   = (round(zper,3),'Standard deviation of photometric zero-point')
	hdr['Starnum'] = (len(zplist_alive),'The number of stars used to obtain zero-point')
	fits.writeto(inim, data, header=hdr, overwrite=True)
	print("FITS Keywords 'Seeing', 'ZP', 'ZPstd', 'Starnum' are updated.")

	# (3) 

	# [1] All stars matched btw. catalog and image.
	allra, alldec = incat['ra'], incat['dec'] 
	allstar       = incat['col1'] # Number
	allcolor      = "green"
	radius        = """ 3" """
	head1 = "# Region file format: DS9 version 4.1\n"
	head2 = """global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n"""
	head3 = "fk5\n"
	g = open(inim[:-5]+'.all.reg','w')
	g.write(head1)
	g.write(head2)
	g.write(head3)
	for n in xrange(len(allra)):
		body1 = "circle("+str(allra[n])+","+str(alldec[n])+","+radius+") # color="+allcolor+" text={"+str(allstar[n])+' '+str(incat['col4'][n])+' '+str(incat['col5'][n])+ "}\n"
		g.write(body1)
	g.close()

	# [2] index selected stars
	refra, refdec = intbl_alive['ra'], intbl_alive['dec']
	starname      = intbl_alive['col1']
	color         = "yellow"
	f = open(inim[:-5]+'.zp.reg','w')
	f.write(head1)
	f.write(head2)
	f.write(head3)
	for n in xrange(len(refra)):
		body2 = "circle("+str(refra[n])+","+str(refdec[n])+","+radius+") # color="+color+" text={"+str(starname[n])+' '+str(intbl_alive['col4'][n])+' '+str(intbl_alive['col5'][n])+ "}\n"	
		f.write(body2)
	f.close()

	### SDSS mag VS Zero Point plot
	plt.figure(1) 
	plt.xlabel('$'+band+'$ [AB mag]')
	plt.ylabel('ZP [AB mag]')
	plt.title(inim)
	if 'fscale' not in inim.split('.') :
		plt.xlim(( 16.5, 21.5))
		plt.ylim(( 27, 29.5))
		plt.plot([0, 30],[zp, zp], linestyle='--', color='black', linewidth=2)
		plt.text(17.5, 29.25, 'ZP ('+str(sigma)+'$\sigma$)= '+str(round(zp,3))+'$\pm$'+str(round(zper,3)))
		plt.text(17.5, 29.1, 'Seeing (median) = '+str(round(seeing_dum,3)))
		plt.text(17.5, 28.95, '# of stars = '+str(len(zplist_alive)))	
	elif 'fscale' in inim.split('.') :
		plt.xlim(( 16.5, 21.5))
		plt.ylim(( 29.5, 32))
		plt.plot([0, 30],[zp, zp], linestyle='--', color='black', linewidth=2)
		plt.text(18.0, 30.8, 'ZP ('+str(sigma)+'$\sigma$)= '+str(round(zp,3))+'$\pm$'+str(round(zper,3)))
		plt.text(18.0, 30.65, 'Seeing (median) = '+str(round(seeing_dum,3)))
		plt.text(18.0, 30.50, '# of stars = '+str(len(zplist_alive)))	
	plt.plot(nomadmag[index4], zplist, linestyle='None', marker='x', color= 'gray')
	plt.plot(intbl_alive[infilter], zplist_alive, linestyle='None', marker='s', color='red', mfc='none', mec='red', mew = 1.0)
	plt.grid(True)
	plt.savefig(inim[:-5]+'.zp.png', dpi=300)
	plt.close()
'''
import glob
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
mao_list = glob.glob('coadd*.R.fscale.fits')
mao_list.sort()
seeing = []
for i in xrange(len(mao_list)): 
	print mao_list[i]
	MAO_fwhmzp_sdss12(mao_list[i], lowermag = 17, uppermag = 22, pixscale = 0.267, gain = 1.47)
	hdr = fits.getheader(mao_list[i])
	seeing.append(hdr['seeing'])
print 'All done.'	

### Seeing histogram
if hist == True :
	plt.figure(2, figsize=(8,6)) 
	plt.xlabel("""FWHM ["]""", fontsize=17)
	plt.ylabel('#', fontsize=17)
	plt.hist(seeing, bins = np.arange(0.8, 2, 0.1),  alpha=0.6, facecolor='blue', histtype = 'barstacked', edgecolor='black', linewidth=1.2)
	#plt.xlim(0.8, 1.7)
	#plt.ylim(0, 8)
	plt.tick_params(labelsize=14)
	plotname = 'mao.fwhm.coadd.hist.sdss.SE.pdf' 
	plt.savefig(plotname)
	plt.close()
	os.system('evince '+plotname+' &')
'''


