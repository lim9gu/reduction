def astrometry(imlist_name, pixscale) :
	"""
	1. Description 
	: Solving WCS coordinates using Astrometry.net software. For better performance in especially B band images, --use-sextractor mode is added. This mode needs SExtractor configuration files. So please posit configuration files for your working directory. cpulimit 300 is also added to prevent too long processing time for bad images.  

	2. Usage
	>>> astrometry('fd*.fits', obs.SAO.pixscale)

	3. History
	2018.03	Created by G.Lim.
	2018.12.18 Edited by G.Lim. SExtractor mode is added.
	2018.12.21 Edited by G.Lim. Define SAO_astrometry function.
	2019.08.07 Edited by G.Lim. Input list edited. Pixscale is required to consider any pixscale of other telescopes.
	"""
	import glob
	import os, sys
	import subprocess
	import numpy as np
	import lgpy.obs as obs
	from lgpy.hdrcheck import wcscenter
	imlist = glob.glob(imlist_name)
	imlist.sort()
	sexconfig = '/data1/code/astrom.config/astrometry.net.sex'
	print('Solving WCS using Astrometry.net...')
	for i in range(len(imlist)):
		inim = imlist[i]
		#radd, decdd = wcscenter(inim)
		radd       = float(input('ra center [deg] = '))
		decdd      = float(input('dec center [deg] = '))
		com='solve-field '+inim+' --cpulimit 300 --overwrite --use-sextractor  --sextractor-config '+sexconfig+' --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_AUTO --sort-ascending --scale-unit arcsecperpix --scale-low '+str(float(pixscale)-0.1)+' --scale-high '+str(float(pixscale)+0.1)+' --no-remove-lines --uniformize 10 --no-plots  --new-fits a'+inim+' --temp-dir . --ra '+str(radd)+' --dec '+str(decdd)+' --radius 0.5' 
		
		print(str(i)+' th of '+str(len(imlist)))
		os.system(com) 
	orinum = subprocess.check_output('ls '+imlist_name+' | wc -l', shell=True)
	resnum = subprocess.check_output('ls a'+imlist_name+' | wc -l', shell=True)
	print("from "+str(orinum[:-1])+" files , "+str(resnum[:-1])+" files are solved.")
	print("All done.")
	os.system('rm tmp*')
	os.system('rm *.wcs *.rdls *.corr *.xyls *.solved *.axy *.match ')
	print('Astrometry process is complete.')

