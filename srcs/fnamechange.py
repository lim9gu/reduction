def fnamechange(imlist_name, ccd='SAO_STX16803') :
	"""
	1. Description 
	: Change file name of WCS solving images (afd*.fits) using naming sequence following: Calib-SAO-CCD_name-OBJECT-UTDATE-UTSTART-FILTER-EXPTIME.fits. Then the code copies changed files to IMSNGgalaxies directory.

	2. Usage
	>>> fnamechange('af*.fits', ccd='SAO_STX16803')

	3. History
	2018.03    Created by G.Lim.
	2018.12.21 Edited by G.Lim. Define SAO_fnamechange function.
    2019.08.07 Edited by G.Lim
	"""
	import glob
	import os, sys
	import subprocess
	import numpy as np 
	from astropy.io import fits
	lists = glob.glob(imlist_name)
	lists.sort()
	for i in range(len(lists)):
		hdr     = fits.getheader(lists[i])
		EXPTIME = int(hdr['exptime'])
		if ccd == 'KL400' :
			FILTER = input('Filter (B,V,R,I) ? : ')
		elif ccd == 'SAO_STX16803' :
			FILTER = hdr['filter']
		OBJECT = hdr['object']
		if (OBJECT == 'M51') | (OBJECT == 'M51a') :
			print('Object name is corrected.')
			OBJECT = 'M51A'
		UTDATE   = hdr['date-obs'][0:10]
		UTSTART  = hdr['date-obs'][11:19]
		newimage = 'Calib-'+ccd+'-'+OBJECT+'-'+str(UTDATE[0:4])+str(UTDATE[5:7])+str(UTDATE[8:10])+'-'+str(UTSTART[0:2])+str(UTSTART[3:5])+str(UTSTART[6:8])+'-'+FILTER+'-'+str(EXPTIME)+'.fits'
		os.system('/usr/bin/cp '+lists[i]+' '+newimage)
		print('Copy '+lists[i]+' to '+newimage+'.')
	print("Basic preprocessing is finished.")