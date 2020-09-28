###### KCT data processing pipeline 
###### Developed by G.Lim (lim9gu@gmail.com)
###### module edited on 2020.03.06
###### run /home/lim9/anaconda3/lib/python3.7/site-packages/lgpy/KCT/STX16803.py
#=================================================================
def fileset(imlist_name, camera = 'STX16803'):
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
    2019. 06. 06 bias, dark, flat, sci images are classified as header  information 'IMAGETYP'. and if there are no calibrations, the code calls recent master images.
                 This function had an issue that cannot call exptime and filter which don't exist. 'hdrcheck' function is added to calculate MJD and to check if object name is entered correctly by comparing CRVAL1, CRVAL2 with IMSNG galaxy catalog. 
    2020. 02. 24 Image information is integrated as table of "cat". Remove all the space in object and filename. 
    2020. 03. 01 process.py is divided into SAO module including STX16803.py and KL4040.py for each CCD and sCMOS.
    2020. 03. 06 Edited for KCT (CDK14) 
    2020. 03. 10 Minor bugs are fixed. G. Lim
    """    
    import glob
    import os, sys
    import numpy as np
    from pyraf import iraf
    from astropy.io import fits
    from astropy.table import Table
    savedir = '/data1/KCT/'
    # Image size cut
    allimage = glob.glob(imlist_name)
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
        newim = im.replace(" ","")
        os.system("mv '"+im+"' "+newim)
        hdr = fits.getheader(newim)
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
        iraf.hedit(newim, fields='object', value=hdr['OBJECT'].replace(" ",""), verify='no')
        
    cat = Table({'allimage' : allimage, 'XBINNING': xbin, 'YBINNING' : ybin,'IMAGETYP' : IMAGETYP, 'EXPTIME' : EXPTIME, 'FILTER' : FILTER, 'OBJECT' : OBJECT}, names = ['allimage', 'XBINNING', 'YBINNING', 'IMAGETYP', 'EXPTIME', 'FILTER', 'OBJECT'])

    bin1 = np.where((cat['XBINNING'] == 1) & (cat['YBINNING'] == 1))[0]
    bin2 = np.where((cat['XBINNING'] == 2) & (cat['YBINNING'] == 2))[0]
    if len(bin2) !=0 :
        print('There are 2x2 binned data...')
        os.system('/usr/bin/mkdir bin2')
        os.system('/usr/bin/mv '+" ".join(cat['allimage'][bin2])+' ./bin2')    
        print('Find 2x2 binned data at bin2 folder.')
    # Name change
    '''
    os.system('/usr/bin/rename zero bias *tion*zero.fit')
    os.system('/usr/bin/rename Bias bias *tion*Bias.fit')
    os.system('/usr/bin/rename Cal cal Cal*bias.fit')
    os.system('/usr/bin/rename D dk *tion*D*.fit')
    #os.system('rename d dk *tion*?d*.fit')
    os.system('/usr/bin/rename dark dk *tion*dark*.fit')
    os.system('/usr/bin/rename Cal cal Cal*dk*.fit')
    '''
    #calibrations = glob.glob('*.fit')
    #calibrations.sort()
    imbin1 = cat[bin1]['allimage']
    bias_idx = np.where(cat['IMAGETYP'] == 'Bias Frame')[0]
    dark_idx = np.where(cat['IMAGETYP'] == 'Dark Frame')[0]

    '''
    os.system('/usr/bin/rename Skyflat skyflat Skyflat*.fit')
    os.system('/usr/bin/rename SKYFLAT skyflat SKYFLAT*.fit')
    os.system('/usr/bin/rename SkyFlat skyflat SkyFlat*.fit')
    os.system('/usr/bin/rename domeflat Domeflat domeflat*.fit')
    os.system('/usr/bin/rename DomeFlat Domeflat DomeFlat*.fit')
    os.system('/usr/bin/rename DOMEFLAT Domeflat domeflat*.fit')
    os.system('/usr/bin/rename t-000 t_000 Domeflat*.fit')
    os.system('/usr/bin/rename flat_ flat- *flat*.fit')
    '''

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
    if  len(bias) == 0 :
        print('No bias today or bias directory is already created!')
        masterbias_loc = savedir+'masterbias/'
        os.system('/usr/bin/ls '+masterbias_loc)
        masterbias_list = glob.glob(masterbias_loc+'*.fits')
        masterbias_list.sort()
        recent_masterbias = masterbias_list[-1]
        if recent_masterbias != [] :
            print('Recent master bias is found.')
            os.system('/usr/bin/cp '+recent_masterbias+' ./')
            print(recent_masterbias+' is copied.')
        elif recent_masterbias  == [] :
            print('No recent bias.')
            pass
    elif len(bias) != 0 :
        os.system('/usr/bin/mkdir bias')
        biasjoin = " ".join(bias)
        os.system('/usr/bin/mv '+biasjoin+' ./bias')

    sciexptime = list(set(cat[sci_idx]['EXPTIME']))
    darkexptime = list(set(cat[dark_idx]['EXPTIME']))
    print('Check today dark exposure...')
    masterdark_loc = savedir+'masterdark/'
    for exp in sciexptime :
        if exp in darkexptime :
            print(str(int(exp)) + ' O' )
        elif exp not in darkexptime :
            print(str(int(exp))+' X' )
            
            os.system('/usr/bin/ls '+masterdark_loc)
            masterdark_list = glob.glob(masterdark_loc+'2*dark'+str(int(exp))+'.fits')
            masterdark_list.sort()
            
            if len(masterdark_list) >= 1 :
                recent_masterdark = masterdark_list[-1]
                print('Recent master dark of '+recent_masterdark+' is found.')
                os.system('/usr/bin/cp '+recent_masterdark+' ./')    
                print(recent_masterdark+' is copied.')
            elif len(masterdark_list) == 0 :
                print('No masterdark. Stop here :(')
                pass
            '''
            if len([recent_masterdark]) == 1 :
                print('Recent master dark of '+recent_masterdark+' is found.')
                os.system('/usr/bin/cp '+recent_masterdark+' ./')    
                print(recent_masterdark+' is copied.')
            elif len([recent_masterdark]) == 0 :
                print('No masterdark. Stop here :(')
            '''
    if len(dark) != 0:
        os.system('/usr/bin/mkdir dark')
        darkjoin = " ".join(dark)
        os.system('/usr/bin/mv '+darkjoin+' ./dark')    
    elif len(dark) == 0 :
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
            masterflat_loc = savedir+'masterflat_'+band+'/'
            os.system('/usr/bin/ls '+masterflat_loc)
            masterflat_list = glob.glob(masterflat_loc+'2*_n'+band+'flat.*.fits')
            masterflat_list.sort()
            
            if len(masterflat_list) >= 1:
                recent_masterflat = masterflat_list[-1]
                print('Recent master flat of '+recent_masterflat+' is found.')
                os.system('/usr/bin/cp '+recent_masterflat+' ./')
            elif len(masterflat_list) == 0:
                print('No masterflat. Stop here :(')
                pass
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
def biascom(imlist_name, camera='STX16803'):
    """
    1. Description 
    : This function makes a master bias image of KCT using Pyraf. Put bias images to 20XX-XX-XX/bias/ directory. Run this code on 20XX-XX-XX directory, then pyraf chdir task enter bias directory and makes process. Output image is zero.fits and it will be copied on upper directory. Due to iraf.chdir task, you should reset python when this code is finished in order to make confusion of current directory between iraf and python! 
    
    2. Usage 
    : Start on 2018-XX-XX directory. Make bias directory which contains each bias frame. Naming of each bias images should be zero*.fit. Then just use biascom()
    
    >>> biascom()

    3. History
    2018.03    Created by G.Lim.
    2018.12.17 Edited by G.Lim. Define SAO_biascom function. 
    2019.02.07 Assign archive of masterbias in each date by G. Lim
    2020.03.01 Modified for KL4040 process
    2020.03.06 Modified for KCT STX16803 process
    """
    import glob
    import os, sys
    import numpy as np
    from pyraf import iraf
    from astropy.io import fits
    savedir = '/data1/KCT/'
    curdir = os.getcwd()
    '''
    yy   = curdir.split('/')[-1].split('-')[0]
    mm   = curdir.split('/')[-1].split('-')[1]
    dd   = curdir.split('/')[-1].split('-')[2]
    curdate = yy+mm+dd
    '''
    curdate = curdir.split('/')[-1]
    iraf.noao()
    iraf.imred()
    iraf.ccdred()
    iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
    iraf.chdir('./bias')
    input_name = 'bias.list'
    output_name = curdate+'_zero.fits'
    #os.system('ls cal*bias.fit > '+input_name)
    calibrations = glob.glob(imlist_name)
    f    = open(input_name,'w+')
    for i in range(len(calibrations)):
        hdr      = fits.getheader(calibrations[i])
        IMAGETYP = hdr['IMAGETYP']
        if IMAGETYP == 'Bias Frame' :    
            f.write(calibrations[i]+'\n')
    f.close()
    print('Zerocombine is running...')
    iraf.zerocombine(input='@'+input_name, output=output_name, combine='median', reject='minmax', process='no', scale='none', ccdtype='' )
    print('Output master '+output_name+' is created.')
    os.system('/usr/bin/cp '+output_name+' ../')
    os.system('mkdir '+savedir+'masterbias')
    os.system('/usr/bin/cp '+output_name+' '+savedir+'masterbias/')
    iraf.chdir('../')
    iraf.dir('.')
#=================================================================
def darkcom(imlist_name,camera='KL4040') :
    """
    1. Description 
    : This function makes a master dark image of KCT using Pyraf. Put dark images to 201X-XX-XX/dark/ directory. Run this code on 201XXXXX directory, then pyraf chdir task enter dark directory and makes process. Output image is darkXXX.fits (XXX is exposure time. This function will classify each exposure of dark frames!) and it will be copied on upper directory. Due to iraf.chdir task, you should reset python when this code is finished in order to make confusion of current directory between iraf and python! 

    2. Usage
    : Start on 2018-XX-XX directory. Make dark directory which contains each dark frame. Naming of each dark image should be dark*.fit. Then just use darkcom().

    >>> darkcom()

    3. History
    2018.03    Created by G.Lim.
    2018.12.20 Edited by G.Lim. Define SAO_darkcom function.
    2019.02.07 Assign archive of masterdark in each date by G. Lim
    2020.03.01 Modified for KL4040 process
    2020.03.06 Modified for KCT STX16803 process
    """
    import glob
    import os, sys
    import numpy as np    
    from pyraf import iraf
    from astropy.io import fits
    savedir = '/data1/KCT/'
    curdir = os.getcwd()
    '''
    yy   = curdir.split('/')[-1].split('-')[0]
    mm   = curdir.split('/')[-1].split('-')[1]
    dd   = curdir.split('/')[-1].split('-')[2]
    curdate = yy+mm+dd
    '''
    curdate = curdir.split('/')[-1]
    iraf.noao()
    iraf.imred()
    iraf.ccdred()
    iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
    iraf.chdir('./dark')

    print('zero subtraction...')
    os.system('ls dark*.fit > dark.list')
    iraf.imarith(operand1='@dark.list', op='-', operand2='../*_zero.fits', result='z@dark.list')
    zdark = glob.glob('zdark*.fit')
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
        input_name  = output_name[:-5]+'.list'
        f=open(input_name,'w+')
        for k in range(len(imlist)) : 
            f.write(imlist[k]+'\n')
        f.close()
        print('Darkcombine is running...')
        iraf.imstat('@'+input_name)
        iraf.darkcombine(input='@'+input_name, output=output_name, combine='median', reject='minmax', process='no', scale='none', ccdtype='' )
    os.system('/usr/bin/cp '+output_name+' ../')
    os.system('mkdir '+savedir+'masterdark')
    os.system('/usr/bin/cp '+output_name+' '+savedir+'masterdark/')
    os.system('/usr/bin/rm d*.list')
    iraf.chdir('../')
    iraf.dir('.')
    print('Output master '+output_name+' is created.')
#=================================================================
def flatcom(camera='STX16803', flattype='sky', dark=False) :
    """
    1. Description 
    : This function makes master-normalised images of KCT using Pyraf. Put flat images to 20XX-XX-XX/skyflat/ directory. Run this code on 20XX-XX-XX directory, then pyraf chdir task enter skyflat, directory and makes process. If dark=True, dark subtraction will perform on flat images. Use this keyword when you think dark subtraction is needed for flat images. If not, only bias subtraction will be perfromed. And then flatcombine and normalizing will be performed. Output image is nflatX.YYY.fits (X is filter and YYY is sky or dome. This function will classify each exposure of frames!) and they will be copied on upper directory. Due to iraf.chdir task, you should reset python when this code is finished in order to make confusion of current directory between iraf and python! 

    2. Usage
    : Start on 20XX-XX-XX directory. Make skyflat or domeflat directory which contains each flat frame. Naming of each flat image should be *flat*.fit. And domeflat naming is Domeflat*.fit. Then just use flatcom(). 
    >>> flatcom('sky') --> Use skyflat
    >>> flatcom('dome') --> Use domeflat
    *Default configuration is skyflat.

    3. History
    2018.03    Created by G. Lim.
    2018.12.20 Edited by G. Lim. Define SAO_flatcom function.
    2018.12.28 Edited by G. Lim. Add bias or dark keyword. Join function is used when performing imarith, combine tasks.
    2019.02.07 Assign archive of masterflat in each date by G. Lim
    2020.03.01 Remove process keyword from STX16803. 
    2020.03.06 Modified for KCT STX16803 process
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
    '''
    yy   = curdir.split('/')[-1].split('-')[0]
    mm   = curdir.split('/')[-1].split('-')[1]
    dd   = curdir.split('/')[-1].split('-')[2]
    curdate = yy+mm+dd
    '''
    curdate = curdir.split('/')[-1]
    savedir = '/data1/KCT/'
    if flattype == 'dome' :
        iraf.chdir('./domeflat')
    elif flattype == 'sky' :
        iraf.chdir('./skyflat')
    flat = glob.glob('*flat*.fit')
    flat.sort()
    input_name = 'flat.list'
    k=0
    f=open(input_name,'w+')
    for k in range(len(flat)) : 
        f.write(flat[k]+'\n')
    f.close()
    print('zero subtraction with '+flattype+'flat images...')
    iraf.imarith(operand1='@'+input_name, op='-', operand2='../*_zero.fits', result='z@'+input_name)
    if dark == True :
        zflat = glob.glob('z*flat*.fit')
        i=0
        allexptime = []
        for i in range(len(zflat)) :
            hdr = fits.getheader(zflat[i])
            allexptime.append(hdr['exptime'])
        expset = set(allexptime)
        exptime = list(sorted(expset))
        i=0
        for i in range(len(exptime)) :
            print('Find images with exptime of '+str(int(exptime[i])))
            imlist = []
            j=0
            for j in range(len(zflat)) :
                hdr = fits.getheader(zflat[j])
                if hdr['exptime'] == exptime[i] :
                    imlist.append(zflat[j])
                else :
                    pass
            print(imlist)
            imlist.sort()
            input_name = 'zflat.list'
            k=0
            f=open(input_name,'w+')
            for k in range(len(imlist)) : 
                f.write(imlist[k]+'\n')
            f.close()
            iraf.imarith(operand1='@'+input_name, op='-', operand2='../*_dark'+str(int(exptime[i]))+'.fits', result='d@'+input_name)
            #iraf.imarith(operand1=darkjoin, op='-', operand2='../dark'+str(int(exptime[i]))+'.fits', result=outdarkjoin)
            print('Dark subtracted flat images are created.')
    # Flat combine
    if dark == True :
        calflat = glob.glob('dz*flat*.fit')
    elif dark == False :
        calflat = glob.glob('z*flat*.fit')
    calflat.sort()
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
        calimlist = []
        for j in range(len(calflat)) :
            hdr = fits.getheader(calflat[j])
            if hdr['filter'] == infilter[i] :
                calimlist.append(calflat[j])
            else :
                pass
        print(calimlist)
        calimlist.sort()
        input_name = str(infilter[i])+'flat.list'
        k=0
        f=open(input_name,'w+')
        for k in range(len(calimlist)) : 
            f.write(calimlist[k]+'\n')
        f.close()

        output_name = input_name[:-5]+'.fits'
        iraf.flatcombine(input='@'+input_name, output=output_name, combine='average', reject='crreject', process='no', scale='mode', ccdtype='', lsigma='3.', hsigma='3.' )

        print(output_name+' is created. Normalizing...')
        data, newhdr = fits.getdata(output_name, header=True)
        x = np.mean(data)
        nimage = data/x
        newflat_name = curdate+'_n'+str(infilter[i])+'flat.'+flattype+'.fits'
        fits.writeto(newflat_name, nimage, header=newhdr, overwrite=True)
        os.system('/usr/bin/cp '+newflat_name+' ../')
        os.system('mkdir '+savedir+'masterflat_'+infilter[i]+'/')
        os.system('/usr/bin/cp '+newflat_name +' '+savedir+'masterflat_'+infilter[i]+'/')
    print('Normalised master flats are created.')
    iraf.imstat(images='*n?flat.'+flattype+'.fits')
    os.system('/usr/bin/rm *.list ?flat.fits')
    iraf.chdir('../')
    iraf.dir('./')
#=================================================================
def objpre(sci_list='*-00*.fit') :
    """
    1. Description
    : This function applies master calibration images to science frames, including bias, dark subtraction, flat fielding. 
 
    2. Usage
    : objpre('*-00*.fit')

    3. History
    2018.03    Created by G. Lim
    2019.02.07 Change name of master calibration frames in each date by G. Lim
    2020.03.01 Bias subtraction is added for KL4040.
    2020.03.06 Modified for KCT STX16803 process    
    """
    import glob
    import os, sys
    import itertools
    import numpy as np
    from pyraf import iraf
    from astropy.io import fits
    from astropy.io import ascii
    curdir = os.getcwd()
    '''
    yy   = curdir.split('/')[-1].split('-')[0]
    mm   = curdir.split('/')[-1].split('-')[1]
    dd   = curdir.split('/')[-1].split('-')[2]
    curdate = yy+mm+dd
    '''
    curdate = curdir.split('/')[-1]
    savedir = '/data1/KCT/'
    iraf.noao()
    iraf.imred()
    iraf.ccdred()
    iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
    #obj_list = '*-00*.fit'
    #obj = glob.glob(obj_list)
    
    # Bias subtraction
    os.system('ls '+sci_list+' > sci.list')
    print('zero subtraction with skyflat images...')
    iraf.imarith(operand1='@sci.list', op='-', operand2='./*_zero.fits', result='z@sci.list')

    obj = glob.glob('z'+sci_list)
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
def astrometry(imlist_name='fd*.fit', scalelow=0.72, scalehigh=0.73, overwrite=False) :
    """
    1. Description 
    : Solving WCS coordinates using Astrometry.net software. For better performance in especially B band images, --use-sextractor mode is added. This mode needs SExtractor configuration files. So please posit configuration files for your working directory. cpulimit 300 is also added to prevent too long processing time for bad images.  

    2. Usage
    >>> astrometry()

    3. History
    2018.03    Created by G.Lim.
    2018.12.18 Edited by G.Lim. SExtractor mode is added.
    2018.12.21 Edited by G.Lim. Define SAO_astrometry function.
    2020.03.01 --backend-config is added to have the system find INDEX files.
    """
    import os,sys
    import glob
    import subprocess
    import numpy as np
    addlist = glob.glob(imlist_name)
    addlist.sort()
    sexconfig = '/data1/code/astrom.config/astrometry.net.sex'
    print('Solving WCS using Astrometry.net...')
    for n in range(len(addlist)):
        if overwrite == False :
            com='/usr/bin/solve-field '+addlist[n]+' --backend-config /usr/local/astrometry/etc/backend.cfg --cpulimit 300 --overwrite --use-sextractor  --sextractor-config '+sexconfig+' --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_AUTO --sort-ascending --scale-unit arcsecperpix --scale-low '+str(scalelow)+' --scale-high '+str(scalehigh)+' --no-remove-lines --uniformize 0 --no-plots  --new-fits a'+addlist[n]+' --temp-dir .\n'
        elif overwrite == True :
            com='/usr/bin/solve-field '+addlist[n]+' --backend-config /usr/local/astrometry/etc/backend.cfg --cpulimit 300 --overwrite --use-sextractor  --sextractor-config '+sexconfig+' --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_AUTO --sort-ascending --scale-unit arcsecperpix --scale-low '+str(scalelow)+' --scale-high '+str(scalehigh)+' --no-remove-lines --uniformize 0 --no-plots  --new-fits '+addlist[n]+' --overwrite --temp-dir .\n'

        print(com)
        print(str(n)+' th of '+str(len(addlist)))
        os.system(com) 
    orinum = subprocess.check_output('ls '+imlist_name+' | wc -l', shell=True)
    resnum = subprocess.check_output('ls a'+imlist_name+' | wc -l', shell=True)
    print("from "+str(orinum[:-1])+" files , "+str(resnum[:-1])+" files are solved.")
    print("All done.")
    os.system('rm tmp*')
    os.system('rm *.wcs *.rdls *.corr *.xyls *.solved *.axy *.match ')
    print('Astrometry process is complete.')
#=================================================================
def hdrcheck(imlist_name='a*.fit', camera='KCT_STX16803'):
    """
    1. Description 
    : When performing observation, header input could be entered wrong. This issue can name the calibrated files inconsistent way for IMSNG survey. To resolve this problem, header ['OBJECT'] would be modified to the name in IMSNG target catalog (alltarget.dat) when the image center is < fov/2 of KCT STX16803. In addition, MJD will be entered, so even if input images are not IMSNG target, they should be processed by this function.

    2. Usage
    >>> hdrcheck()

    3. History
    2018.      Firstly made.
    2020.03.01 Edited for KL4040.
    2020.03.06 Modified for KCT STX16803 process    

    """
    import os
    import sys
    import glob
    import astropy.units as u
    from astropy.time import Time
    from astropy.io import ascii
    from astropy.io import fits
    from lgpy.hdrcheck import wcscenter
    from astropy.coordinates import SkyCoord

    KCT_fov = 49.4 # arcmin
    all_catname = '/data1/code/alltarget.dat'
    all_cat     = ascii.read(all_catname) 
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
        inim            = imlist[i]
        print(inim)
        data, hdr       = fits.getdata(inim, header=True)
        CRVAL1, CRVAL2  = wcscenter(inim)

        #camera = inim.split('-')[1]
        mjd0 = 2400000.5
        if camera == 'MAO_SNUCAM' :
            t = Time(hdr['utdate']+'T'+hdr['utstart'], format='isot', scale='utc')
            jd = t.jd
            mjd = t.mjd
            hdr['JD'] = round(jd,5)
        elif camera == 'KCT_STX16803' :
            t = Time(hdr['DATE-OBS'], format='isot', scale='utc')
            jd = t.jd
            mjd = t.mjd
        else: 
            jd = hdr['jd']
            mjd  = jd - mjd0
        hdr['MJD']  = round(mjd,5)
        coo_target      = SkyCoord(CRVAL1, CRVAL2, unit=(u.deg, u.deg))
        indx, d2d, d3d  = coo_target.match_to_catalog_sky(coo_all)

        if d2d.arcmin > KCT_fov/2. :
            print('Coordinates of the image are not in IMSNG catalog. No matching. Maybe you obtained wrong field. OR Non-IMSNG target.' )
            fits.writeto(inim, data, header=hdr, overwrite=True)
            print('Only MJD is entered in image header.')
            pass
        elif d2d.arcmin < KCT_fov/2. :
                obj  = all_cat[indx]['obj']
                print('======================================')
                print(obj+ ' is matched.')
                print(str(round(d2d.arcmin[0],3))+ ' arcmin apart')
                print('======================================')
                hdr['object']   = obj
                fits.writeto(inim, data, header=hdr, overwrite=True)
                
    print('Header info inspection is finished.')
#=================================================================
def fnamechange(imlist_name,camera='KCT_STX16803') :
    """
    1. Description 
    : Change file name of WCS solving images (afd*.fits) using naming sequence following: Calib-KCT_CCD_name-OBJECT-UTDATE-UTSTART-FILTER-EXPTIME.fits. Then the code copies changed files to IMSNGgalaxies directory.

    2. Usage
    >>> fnamechange()

    3. History
    2018.03    Created by G.Lim.
    2018.12.21 Edited by G.Lim. Define SAO_fnamechange function.
    2020.03.01 Edited for KL4040.
    2020.03.06 Modified for KCT STX16803 process    
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
        elif camera == 'SAO_KL4040' :
            OBJECT = hdr['object']
           
            FILTER = hdr['filter']
        elif camera == 'SAO_STX16803' :
            OBJECT = hdr['object']
            FILTER = hdr['filter']
        elif camera == 'CCA250' :
            OBJECT = hdr['object']
            FILTER = hdr['filter']
        elif camera == 'KCT_STX16803' :
            #OBJECT = hdr['object']
            OBJECT = hdr['OBJECT'].replace(" ","")
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
    print("Basic preprocessing of KCT is finished.")
#=================================================================
def filemove( password, camera='KCT_STX16803') :
    """
    1. Description 
    : This code distributes all the calibrated images of KCT data to IMSNGgalaxies/KCT directory, based on C. Choi's code.
    
    2. Usage
    : Run this code on '/data3/IMSNG/IMSNGgalaxies' location.
    >>> filemove() 

    3. History
    2018.12    Created by G.Lim 
    2018.01.24 Docstring is added by G. Lim
    2020.03.01 For my desktop, move calibrated images to the server using scp.
    2020.03.06 Modified for KCT STX16803 process    
    """
    import os
    import numpy as np
    import astropy.io.fits as fits
    import astropy.io.ascii as ascii
    loc = '/data3/IMSNG/IMSNGgalaxies/'
    os.system('ls Cal*.fits -l > currentdir.list')
    curdirlist = ascii.read('currentdir.list',data_start=0)
    dirpar     = curdirlist['col2']
    name       = curdirlist['col9']
    #dirnames   = name[np.where(dirpar >  1)]
    filenames  = name[np.where(dirpar == 1)]
    calframes  = []
    for i in filenames : 
        if (i[:9] =='Calib-KCT') & (i[-4:] == 'fits') : calframes.append(i)
    print(str(len(calframes))+' files exist.')
    for n in calframes:
        galname  = n.split('-')[2]
        frontname2 = list(galname)[0]+ list(galname)[1]
        frontname3 = list(galname)[0]+ list(galname)[1]+ list(galname)[2]
        if (frontname3  == 'ESO') | (frontname3  == 'KUG') | (frontname3  == 'PKS'):
            galname  = n.split('-')[2]+'-'+n.split('-')[3]
        elif (frontname2  == 'AM') :
            galname  = n.split('-')[2]+'-'+n.split('-')[3]
        else :
            pass
        print(galname)
        if galname == 'M51a' :
            galname = 'M51A'
        hdr = fits.getheader(n)
        if camera == 'STX16803' :
            infilter = hdr['filter']
        elif camera == 'KL4040' :
            infilter = hdr['filter']
        elif camera == 'KCT_STX16803' :
            infilter = hdr['filter']
        elif camera == 'Kepler' :
            infilter = 'R'
        if password == '':
            input('Enter the password for your server account : ')
        else :
            makedir    = 'sshpass -p "'+password+'" ssh lim9@qso.snu.ac.kr /usr/bin/mkdir -p '+loc+galname
            os.system(makedir)
            makeobsdir = 'sshpass -p "'+password+'" ssh lim9@qso.snu.ac.kr /usr/bin/mkdir -p  '+loc+galname+'/KCT/' # If IMSNG target
            os.system(makeobsdir)
            makecamdir = 'sshpass -p "'+password+'" ssh lim9@qso.snu.ac.kr /usr/bin/mkdir -p  '+loc+galname+'/KCT/'+camera
            os.system(makecamdir)
            makefildir = 'sshpass -p "'+password+'" ssh lim9@qso.snu.ac.kr /usr/bin/mkdir -p  '+loc+galname+'/KCT/'+camera+'/'+infilter+'/'
            os.system(makefildir)
            mvcommand='sshpass -p "'+password+'" scp '+ n +' lim9@qso.snu.ac.kr:'+loc+galname+'/KCT/'+camera+'/'+infilter+'/'
            mvsubtract = 'sshpass -p "'+password+'" scp hd'+n+' hc'+n+' lim9@qso.snu.ac.kr:'+loc+galname+'/KCT/'+camera+'/'+infilter+'/'
            os.system(mvcommand)
            os.system(mvsubtract)
            chmodcom = 'sshpass -p "'+password+'" ssh lim9@qso.snu.ac.kr /usr/bin/chmod 777 -R '+loc+galname+'/KCT/'
            os.system(chmodcom)
    print('Files were transferred to the server.')
#=========================================================================
def main(calzp = False):
    '''
    1. Description
    Run basic reduction for KL4040.
    Run this function on 20XX-XX-XX directory.

    2. Usage 
    >>> main()

    3. History
    2020.03.01 Created by G.Lim
    2020.03.06 Edited for KCT STX16803, need flux calibration function for SDSS filters. 
    '''
    import os, sys
    import glob
    from astropy.io import fits
    from lgpy.KCT import STX16803
    from lgpy import photcom
    from lgpy import run_as
    print('Basic reduction starts...')
    os.system('rename fts fit *.fts')
    STX16803.fileset('*-00*.fit',camera='STX16803')
    mb = glob.glob('2*zero.fits')
    if len(mb) == 1 :
        print('I have '+mb[0]+'.')
        pass
    elif len(mb) == 0 :
        STX16803.biascom('zero*.fit',camera='STX16803')
    d = glob.glob('./dark/*.fit')
    if len(d) != 0 :
        STX16803.darkcom(camera='STX16803')
    elif len(d) == 0 :
        md = glob.glob('2*dark*.fits')
        print('I have')
        print(md)
    f = glob.glob('*flat*.fit') 
    if len(f) != 0 :
        STX16803.flatcom(flattype='sky',  camera='STX16803')
    elif len(f) == 0 :
        mf = glob.glob('2*n*.fits')
        print('I have')
        print(mf)      
    STX16803.objpre(sci_list='*-00*.fit')
    run_as.run_as('fdz*.fit', 1.27, 0.724, 9.0, objlim=10 )
    STX16803.astrometry(imlist_name='cfd*.fit')
    STX16803.hdrcheck(imlist_name='a*.fit')
    STX16803.fnamechange('a*.fit',camera='KCT_STX16803')
    print('Basic reduction is finished. ')
    if calzp == True :
        print('calzp = True. Photometry starts...')
        calim = glob.glob('Cal*.fits')
        calim.sort()
        for i in range(len(calim)):
            inim = calim[i]
            band = fits.getheader(inim)['FILTER']
            obj  = fits.getheader(inim)['OBJECT']
            if band == 'u' :
                photcom.run_phot(inim, 'PS1', 'KCT_STX16803', inmagkey='MAG_APER_6', inmagerkey='MAGERR_APER_6', refmagkey='u', refmagerkey='uerr', sigma_n = 5, skylim=True)
            elif band == 'g' :
                if obj == 'ESO182-G010' :
                    photcom.run_phot(inim, 'APASS', 'KCT_STX16803', inmagkey='MAG_APER_6', inmagerkey='MAGERR_APER_6', refmagkey='g', refmagerkey='gerr', sigma_n = 5, skylim=True)
                else :
                    photcom.run_phot(inim, 'PS1', 'KCT_STX16803', inmagkey='MAG_APER_6', inmagerkey='MAGERR_APER_6', refmagkey='g', refmagerkey='gerr', sigma_n = 5, skylim=True)                
            elif band == 'r' :
                if obj == 'ESO182-G010' :
                    photcom.run_phot(inim, 'APASS', 'KCT_STX16803', inmagkey='MAG_APER_6', inmagerkey='MAGERR_APER_6', refmagkey='r', refmagerkey='rerr', sigma_n = 5, skylim=True)  
                else :
                    photcom.run_phot(inim, 'PS1', 'KCT_STX16803', inmagkey='MAG_APER_6', inmagerkey='MAGERR_APER_6', refmagkey='r', refmagerkey='rerr', sigma_n = 5, skylim=True)        
            elif band == 'i' :
                if obj == 'ESO182-G010' :
                    photcom.run_phot(inim, 'APASS', 'KCT_STX16803', inmagkey='MAG_APER_6', inmagerkey='MAGERR_APER_6', refmagkey='i', refmagerkey='ierr', sigma_n = 5, skylim=True)  
                else :
                    photcom.run_phot(inim, 'PS1', 'KCT_STX16803', inmagkey='MAG_APER_6', inmagerkey='MAGERR_APER_6', refmagkey='i', refmagerkey='ierr', sigma_n = 5, skylim=True)
        os.system('gethead Cal*.fits zp zper exptime fwhm_SE lim_SE')
    elif calzp == False :
        print('calzp = False. No photometry.')
    STX16803.filemove('kangsun710', camera='KCT_STX16803')
    print('ZP, FWHM estimation is complete. Subtraction is TBD.')


