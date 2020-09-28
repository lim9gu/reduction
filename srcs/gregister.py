def run_gregister(imlist_name, refim, input_list=False):
	"""
	1. Description 
	: Specify reference image into refim. images to align are imlist_name.
	  Compare pixscale of each input sci image and ref image.
		- Pixscale_sci >= Pixscale_ref : Bin reference image.
		- Pixscale_sci <  Pixscale_ref : Bin science image. ZP should be calculated again.
	2. Usage
	>>> run_gregister('Cal*0.fits', 'Ref.fits', input_list=False)
	>>> run_gregister(np.array(['sci1.fits, sci2.fits, sci3.fits']), 'Ref.fits', input_list=True)

	3. History
	2019.06.13  Created by G. Lim
	2019.08.07  Edited by G. Lim
	"""
	import glob
	import alipy
	import os,sys
	import numpy as np
	from astropy.io import fits
	import lgpy.obs as obs

	def gregister(images):
		import alipy
		id     = alipy.ident.run(ref_image, images, visu=False)
		print("%20s : %20s, flux ratio %.2f" % (id[0].ukn.name, id[0].trans, id[0].medfluxratio))
		alipy.align.irafalign(id[0].ukn.filepath, id[0].uknmatchstars, id[0].refmatchstars, shape=outputshape, makepng=False)	
	if input_list == False :
		image_list_filter = glob.glob(imlist_name)	
	elif input_list == True :
		image_list_filter = list(imlist_name)
	image_list_filter.sort()
	obs_sci    = image_list_filter[0].split('-')[1]

	ref_image = refim
	obs_ref   = ref_image.split('-')[1]
	if obs_ref in ['SDSS', 'PS1'] :
		obs_ref = input('Enter reference obs :')

	if obs.pixscale(obs=obs_sci) < obs.pixscale(obs=obs_ref) :
		if len(image_list_filter) > 1 : 
			print('Pixscale_original < Pixscale_reference')
			images_to_align    = image_list_filter
			ref_image          = refim
			outputshape        = alipy.align.shape(ref_image)
			images_to_align_1  = []
			n = 0
			for n in range(len(images_to_align)) :
				images        = images_to_align[n:n+1][0]
				print('\n',images,'\n')
				gregister([images])
				images_to_align_1.append(images)
				os.system('mv alipy_out/'+images[:-5]+'_gregister.fits '+'alipy_out/gRef'+images)
		else :
			print(str(image_list_filter[0])+' is the only element. Pass.')
			pass
	elif obs.pixscale(obs=obs_sci) >= obs.pixscale(obs=obs_ref) :
		if len(image_list_filter) > 0 : 
			for i in range(len(image_list_filter)):
				print('Pixscale_original >= Pixscale_reference')
				inim    = image_list_filter[i]
				part    = inim.split('-')
				part[0] = 'Ref'
				newref  = '-'.join(part)
				os.system('cp '+refim+' '+newref)
				images_to_align = newref
				ref_image = inim
				print('ref : '+inim) 
				outputshape       = alipy.align.shape(ref_image)
				images_to_align_1 = []
				images        = images_to_align
				gregister([images])
				images_to_align_1.append(images)
				os.system('mv alipy_out/'+newref[:-5]+'_gregister.fits '+'alipy_out/g'+newref)
				os.system('rm '+newref)
		else :
			print(str(image_list_filter[0])+' is the only element. Pass.')
			pass
	os.system('mv alipy_out/gRef*.fits .')
	os.system('rm -r alipy_out')
	print('Done. \a')