def FluxScale(inim_list, zp0):
	"""
	Multiply constant to match zero point of images to combine taken in different sky, weather, exposure time.

	>>> FluxScale(imlist_name, refim=0)

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