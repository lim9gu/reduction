###### KHAO data processing pipeline
###### Developed by G.Lim (lim9gu@gmail.com)
###### module edited on 2019.08.07
###### run /home/lim9/anaconda3/lib/python3.7/site-packages/lgpy/KHAO/process.py
#=================================================================
def fileset(imlist_name='*_00*.fit'):
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
    2019. 06. 06 bias, dark, flat, sci images are classified as header information 'IMAGETYP'. and if there are no calibrations, the code calls recent master images. This function had an issue that cannot call exptime and filter which don't exist. 'hdrcheck' function is added to calculate MJD and to check if object name is entered correctly by comparing CRVAL1, CRVAL2 with IMSNG galaxy catalog. 
    """    
    import glob
    import os, sys
    import numpy as np
    from astropy.io import fits
    from astropy.table import Table
    # Image size cut
    allimage = glob.glob(imlist_name)
    allimage.sort()
    curdir = os.getcwd()
    curdate = curdir.split('/')[-1]
    # Header sorting
    image, xbin, ybin, imtype, band, exptime = [], [], [], [], [], []
    for im in allimage :
        hdr        = fits.getheader(im)
        xbin_dum   = hdr['XBINNING']
        ybin_dum   = hdr['YBINNING']
        imtype_dum = hdr['IMAGETYP']
        band_dum   = hdr['FILTER']
        exptime_dum = hdr['exptime']
        
        xbin.append(xbin_dum)
        ybin.append(ybin_dum)  
        imtype.append(imtype_dum)
        image.append(im)
        band.append(band_dum)
        exptime.append(exptime_dum)

    cat = Table(  {'IMAGE' : image, 'XBIN': xbin, 'YBIN' : ybin, 'IMTYPE' : imtype, 'BAND' : band, 'EXPTIME' : exptime}, names = ['IMAGE', 'XBIN', 'YBIN', 'IMTYPE', 'BAND', 'EXPTIME'])

    # Binning classification
    bin2 = np.where((cat['XBIN'] == 2) & (cat['YBIN'] == 2))[0]
    bin3 = np.where((cat['XBIN'] == 3) & (cat['YBIN'] == 3))[0]
    if len(bin3) !=0 :
        print('There are 3x3 binned data...')
        os.system('/usr/bin/mkdir bin3')
        os.system('/usr/bin/mv '+" ".join(cat['IMAGE'][bin3])+' ./bin3')
        print('Find 3x3 binned data at bin3 folder.')

    # IMAGE CLASSIFICATION
    imbin2 = cat[bin2]['IMAGE']

    bias   = cat[bin2]['IMAGE'][np.where(cat[bin2]['IMTYPE'] == 'BIAS')]
    input_bias = 'bias.list'
    f    = open(input_bias,'w+')
    for i in range(len(bias)):    
        f.write(bias[i]+'\n')
    f.close()

    dark   = cat[bin2]['IMAGE'][np.where(cat[bin2]['IMTYPE'] == 'DARK')]
    input_dark = 'dark.list'
    f    = open(input_dark,'w+')
    for i in range(len(dark)):    
        f.write(dark[i]+'\n')
    f.close()   

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

        input_dark = 'dark'+str(int(exptime[i]))+'.list'
        output_dark = curdate+'_dark'+str(int(exptime[i]))+'.fits'
        #input_name = output_name[:-5]+'.list'
        f=open(input_dark,'w+')
        for k in range(len(imlist)) : 
            f.write(imlist[k]+'\n')
        f.close()

    flat   = cat[bin2]['IMAGE'][np.where(cat[bin2]['IMTYPE'] == 'FLAT')]
    input_flat = 'flat.list'
    f    = open(input_flat,'w+')
    for i in range(len(flat)):    
        f.write(flat[i]+'\n')
    f.close()   

    light  = cat[bin2]['IMAGE'][np.where(cat[bin2]['IMTYPE'] == 'LIGHT')]
    input_light = 'light.list'
    f    = open(input_light,'w+')
    for i in range(len(light)):    
        f.write(light[i]+'\n')
    f.close()   

    print('All images are classified. bias, dark, flat, light are returned arrays')
    return bias, dark, flat, light
    '''
    # IF NO BIAS TODAY,
    if  bias == [] :
        print('No bias today or bias directory is already created!')
        masterbias_loc = '/data3/SAO1m/red/'+camera+'/masterbias/'
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

    # Dark exp time classification
    masterdark_loc = '/data3/SAO1m/red/'+camera+'/masterdark/'
    os.system('/usr/bin/ls '+masterdark_loc)

    # Exptime of sci images
    sciexptime = []
    for i in range(len(sci)) :
        scihdr = fits.getheader(sci[i])
        sciexptime.append(scihdr['exptime'])
        sciexpset = set(sciexptime)
        sciexptime = list(sorted(sciexpset))    

    for exp in sciexptime:
            darkexp = glob.glob('cal*dk'+str(int(exp))+'.fit')
            if len(darkexp) != 0:
                os.system('/usr/bin/mkdir dark')
                darkjoin = " ".join(darkexp)
                os.system('/usr/bin/mv '+darkjoin+' ./dark')    
            elif len(darkexp) == 0 :
                dark_dum_name = '2*_dark'+str(int(exp))+'.fits'
                print('Find master dark of '+dark_dum_name+'...')
                dark_dum = glob.glob(dark_dum_name)
                if len(dark_dum) == 1 :
                    print('I found '+dark_dum_name+'!')
                    print('Work with this masterdark. STOP.')
                elif len(dark_dum) == 0 :
                    print('No dark images of '+str(int(exp)))
                    masterdark_list = glob.glob(masterdark_loc+'2*dark'+str(int(exp))+'.fits')
                    masterdark_list.sort()
                    print(masterdark_list)
                    recent_masterdark = masterdark_list[-1]
                    if len([recent_masterdark]) == 1 :
                        print('Recent master dark of '+recent_masterdark+' is found.')
                        os.system('/usr/bin/cp '+recent_masterdark+' ')
                        print(recent_masterdark+' is copied.')
                    elif len([recent_masterdark]) == 0 :
                        print('No masterdark. Stop here :(')    

    # Flat classification
    os.system('/usr/bin/rename Skyflat skyflat Skyflat*.fit')
    os.system('/usr/bin/rename SKYFLAT skyflat SKYFLAT*.fit')
    os.system('/usr/bin/rename SkyFlat skyflat SkyFlat*.fit')
    os.system('/usr/bin/rename domeflat Domeflat domeflat*.fit')
    os.system('/usr/bin/rename DomeFlat Domeflat DomeFlat*.fit')
    os.system('/usr/bin/rename DOMEFLAT Domeflat domeflat*.fit')
    os.system('/usr/bin/rename t-000 t_000 Domeflat*.fit')
    os.system('/usr/bin/rename flat_ flat- *flat*.fit')

    # Filter of science images

    scifilter = []
    i=0
    for i in range(len(sci)):
        scihdr = fits.getheader(sci[i])
        scifilter.append(scihdr['filter'])
    scifilterset = set(scifilter)
    sciinfilter = list(sorted(scifilterset))
    print(sciinfilter, 'filter set today.')

    for band in sciinfilter :
        masterflat_loc = '/data3/SAO1m/red/'+camera+'/masterflat_'+band+'/'
        os.system('/usr/bin/ls '+masterflat_loc)
        flatband = glob.glob('*flat*'+band+'.fit')
        if len(flatband) != 0 :
            if 'skyflat' in flatband[0].split('-') :
                skyflat = flatband
                os.system('/usr/bin/mkdir skyflat')
                skyflatjoin = " ".join(skyflat)
                os.system('/usr/bin/mv '+skyflatjoin+' ./skyflat')
            elif 'Domeflat' in flatband[0].split('-') :
                domeflat = flatband
                os.system('/usr/bin/mkdir domeflat')
                domeflatjoin = " ".join(domeflat)
                os.system('/usr/bin/mv '+domeflatjoin+' ./domeflat')            
        elif len(flatband) == 0 :
            print('No flat frames today.')
            masterflat_list = glob.glob(masterflat_loc+'2*_n'+band+'flat.*.fits')
            masterflat_list.sort()
            print(masterflat_list)
            recent_masterflat = masterflat_list[-1]
            if len([recent_masterflat]) == 1:
                print('Recent master flat of '+recent_masterflat+' is found.')
                os.system('/usr/bin/cp '+recent_masterflat+' ./')
            elif len([recent_masterflat]) == 0:
                print('No masterflat. Stop here :(')
    print('File setting is done.')
    '''
#=================================================================
def biascom(input_bias='bias.list'):
    """
    1. Description 
    : This function makes a master bias image of SAO 1-m using Pyraf. Put bias images to 201XXXXX/bias/ directory. Run this code on 201XXXXX directory, then pyraf chdir task enter bias directory and makes process. Output image is zero.fits and it will be copied on upper directory. Due to iraf.chdir task, you should reset python when this code is finished in order to make confusion of current directory between iraf and python! 
    
    2. Usage 
    : Start on 2018XXXX directory. Make bias directory which contains each bias frame. Naming of each bias images should be cal*bias.fit. Then just use SAO_biascom()
    
    >>> SAO_biascom()

    3. History
    2018.03    Created by G.Lim.
    2018.12.17 Edited by G.Lim. Define SAO_biascom function. 
    2019.02.07 Assign archive of masterbias in each date by G. Lim
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
    #input_name = 'bias.list'
    output_name = curdate+'_zero.fits'
    #calibrations = input_bias
    #f    = open(input_name,'w+')
    #for i in range(len(calibrations)):    
    #    f.write(calibrations[i]+'\n')
    #f.close()
    print('Zerocombine is running...')
    iraf.imstat(images='@'+input_bias)
    iraf.zerocombine(input='@'+input_bias, output=output_name, combine='median', reject='minmax', process='no', scale='none', ccdtype='' )
    print('Output master '+output_name+' is created.')
    os.system('/usr/bin/cp '+output_name+' /data1/KHAO/MDFTS/red/masterbias/')
    iraf.dir('.')
#=================================================================
def biassub(input_image, input_masterbias) :
    import os, sys
    from pyraf import iraf
    iraf.imarith(operand1='@'+input_image, op='-', operand2=input_masterbias, result='z@'+input_image)
    print('Masterbias is subtracted from the images.')
#=================================================================
def darkcom(input_dark):
    import glob
    import os, sys
    import numpy as np    
    from pyraf import iraf
    from astropy.io import fits
    from astropy.io import ascii
    curdir = os.getcwd()
    curdate = curdir.split('/')[-1]
    #dark = ", ".join(input_dark)
    iraf.noao()
    iraf.imred()
    iraf.ccdred()
    iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
    #iraf.chdir('./dark')
    #dark = glob.glob('cal*dk*.fit')
    #dark = ascii.read('dark.list', guess=True, data_start=0)
    '''
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
        input_name = 'dark'+str(int(exptime[i]))+'.list'
        output_name = curdate+'_dark'+str(int(exptime[i]))+'.fits'
        #input_name = output_name[:-5]+'.list'
        f=open(input_name,'w+')
        for k in range(len(imlist)) : 
            f.write(imlist[k]+'\n')
        f.close()
    '''
    #output_name = curdate+'_dark'+str(int(exptime[i]))+'.fits'
    output_name = curdate+'_'+input_dark[:-5]+'.fits'
    print('Darkcombine is running...')
    iraf.imstat(images='z@'+input_dark)
    iraf.darkcombine(input='z@'+input_dark, output=output_name, combine='median', reject='minmax', process='no', scale='none', ccdtype='' )
    #os.system('/usr/bin/cp '+output_name+' ../')
    os.system('/usr/bin/cp '+output_name+' /data1/KHAO/MDFTS/red/masterdark/')
    #os.system('/usr/bin/rm d*.list')
    #iraf.chdir('../')
    iraf.dir('.')
    print('Output master '+output_name+' is created.') 
#=================================================================
def darksub(input_image, input_masterdark):
    import os, sys
    from pyraf import iraf
    iraf.imarith(operand1='z@'+input_image, op='-', operand2=input_masterdark, result='dz@'+input_image)
    print('Masterdark is subtracted from the images.')   
#=================================================================
def flatcom(input_flat='zFLAT*.fit'):
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
    # Filter classification
    calflat = glob.glob(input_flat)
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
        output_name = input_name[:-5]+'.fits'
        #iraf.flatcombine(input='@'+input_name, output=output_name, combine='average', reject='crreject', process='no', scale='mode', ccdtype='', lsigma='3.', hsigma='3.' )
        iraf.flatcombine(input='@'+input_name, output=output_name, combine='median', reject='minmax', process='no', scale='mode', ccdtype='')
        print(output_name+' is created. Normalizing...')
        data, newhdr = fits.getdata(output_name, header=True)
        x = np.mean(data)
        nimage = data/x
        newflat_name = curdate+'_n'+str(infilter[i])+'flat.fits'
        fits.writeto(newflat_name, nimage, header=newhdr, overwrite=True)
        #os.system('/usr/bin/cp '+newflat_name+' ../')
        os.system('/usr/bin/cp '+newflat_name+' /data1/KHAO/MDFTS/red/masterflat_'+infilter[i]+'/')
    print('Normalised master flats are created.')
    iraf.imstat(images='*n?flat.fits')
#=================================================================
'''
def flatdiv(input_image, input_masterflat):
    import os, sys
    from pyraf import iraf
    iraf.imarith(operand1='dz@'+input_image, op='/', operand2=input_masterflat, result='fdz@'+input_image)
    print('Masterflat is divided from the images.')      
'''
#=================================================================
def objpre(light):
    import glob
    import os, sys
    import itertools
    import numpy as np
    from pyraf import iraf
    from astropy.io import fits
    from astropy.io import ascii

    curdir = os.getcwd()
    curdate = curdir.split('/')[-1]
    biassub(input_image='light.list', input_masterbias = curdate+'_zero.fits')
    #zobj = ['z'+x for x in light]
    #zobj.sort()
    obj = light
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
        darksub(input_image= 'obj'+str(int(exptime[i]))+'.list', input_masterdark = curdate+'_dark'+str(int(exptime[i]))+'.fits')
    # Flat fielding
    #dobj = glob.glob('d'+obj_list)
    dobj = ['dz'+x for x in obj]
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
        print('Performing flat fielding...')
        nflats = glob.glob('2*n'+str(infilter[i])+'flat.fits')[0]
        flattype = nflats.split('.')[1]
        print(nflats)
        #nflats = 'n'+str(infilter[i])+'flat.'+flattype+'.fits'
        iraf.imarith(operand1='@obj'+infilter[i]+'.list', op='/', operand2=nflats, result='f@obj'+infilter[i]+'.list')
        #iraf.imarith(operand1=dimjoin, op='/', operand2=nflats, result=doutimjoin)
    print('Flat fielding is finished. Check the images.')    
#=================================================================
def astrometry(imlist_name='fdz*.fit'):
    import os,sys
    import glob
    import subprocess
    import numpy as np
    addlist = glob.glob(imlist_name)
    addlist.sort()
    sexconfig = '/data1/code/astrom.config/astrometry.net.sex'
    print('Solving WCS using Astrometry.net...')
    for n in range(len(addlist)):
        com='solve-field '+addlist[n]+' --cpulimit 300 --overwrite --use-sextractor  --sextractor-config '+sexconfig+' --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_AUTO --sort-ascending --scale-unit arcsecperpix --scale-low 0.7 --scale-high 0.8 --no-remove-lines --uniformize 0 --no-plots  --new-fits a'+addlist[n]+' --temp-dir .\n'
        #com='/usr/bin/solve-field '+addlist[n]+' --cpulimit 300 --overwrite --use-sextractor --scale-unit arcsecperpix --scale-low 0.7 --scale-high 0.8 --no-remove-lines --uniformize 0 --no-plots  --new-fits a'+addlist[n]+' --temp-dir .\n'        
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
def fnamechange() :
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
    import subprocess
    import numpy as np 
    from astropy.io import fits
    lists = glob.glob('afd*.fit')
    lists.sort()
    for i in range(len(lists)):
        hdr = fits.getheader(lists[i])
        EXPTIME = int(hdr['exptime'])
        OBJECT = hdr['object']
        if (OBJECT == 'M51') | (OBJECT == 'M51a') :
            print('Object name is corrected.')
            OBJECT = 'M51A'
        UTDATE = hdr['date-obs'][0:10]
        UTSTART = hdr['date-obs'][11:19]
        FILTER  = hdr['FILTER']
        newimage = 'Calib-KHAO_MDFTS-'+OBJECT+'-'+str(UTDATE[0:4])+str(UTDATE[5:7])+str(UTDATE[8:10])+'-'+str(UTSTART[0:2])+str(UTSTART[3:5])+str(UTSTART[6:8])+'-'+FILTER+'-'+str(EXPTIME)+'.fits'
        os.system('/usr/bin/cp '+lists[i]+' '+newimage)
        print('Copy '+lists[i]+' to '+newimage+'.')
    #os.system('cp Cal*.fits /data3/IMSNG/IMSNGgalaxies/')
    print("Basic preprocessing of KHAO_MDFTS is finished.")
#=================================================================
def main():
    import os
    import glob
    from lgpy import hdrcheck
    curdir = os.getcwd()
    curdate = curdir.split('/')[-1]
    bias, dark, flat, light = fileset(imlist_name='*_00*.fit')
    biascom(input_bias='bias.list')
    biassub(input_image='dark.list', input_masterbias = curdate+'_zero.fits')
    darklist = glob.glob('dark*.list')
    darklist.sort()
    darklist.remove(darklist[0])
    i=0
    for i in range(len(darklist)) :
        darkcom(input_dark=darklist[i])
    objpre(light)
    astrometry(imlist_name='fdz*.fit')
    hdrcheck(imlist_name='afdz*.fit', ccd = 'KHAO_MDFTS', fov=16., racen='', deccen='')
    fnamechange()