def apall(imlist_name, reference=''):
    """
    Extract 1d spectrum by tracing column of 2d data.
    Interactive mode is turned on by default.

    On the interactive mode;
    (1) Aperture selection
    :w -> e1, e2 : Zoom in
    :w -> a      : Zoom out
    :m           : Select aperture
    :l -> Click left of the aper. : Set lower limit of the aperture.
    :u -> Click right of the aper. : Set upper limit of the aperture.

    (2) Background selection
    :z  : Delete previously-selected background region.
    :b  : Enter the background mode
    :s1, s2 (on the left part of bkg) : Set the region of left bkg.
    :s3, s4 (on the right part of bkg) : Set the region of right bkg.
    And then enter q --> Yes, Yes

    (3) Fitting mode
    :d : Exclude a point
    :function : Change the fitting function. ("legendre" by default) 
    :order=15 : Change the function order into 15.
    :niterate : Cliping iteration (3 times by default)
    :f        : Go fitting
    :q        : Confirm
    And the enter q --> Yes, Yes

    Output image will be saved as "XXX.ms.fits"
    reference='' : Not using specified aperture used.
    reference='??.ms' : Use the same aperture of ??.ms. Do copy of database of the reference into the same location of input image***
    """
    import glob
    import os, sys
    from pyraf import iraf
    iraf.noao()
    iraf.twodspec()
    iraf.apextract()
    imlist = glob.glob(imlist_name)
    imlist.sort()
    for i in range(len(imlist)):
        inim = imlist[i]
        print('1d extraction for '+inim+'...')
        if reference == '' :
            iraf.apall(input=inim, apertures='', interactive='yes', find='no', t_function="legendre", t_order=15, t_niterate=3, t_low_reject = 3, t_high_reject = 3. )
        elif reference != '' :
            iraf.apall(input=inim, apertures='', reference=reference, interactive='no', find='no', recenter='no', trace='no',t_function="legendre", t_order=15, t_niterate=3, t_low_reject = 3, t_high_reject = 3. )            
    print('Done.')

def rotate(imlist_name, rotation):
    """
    rotate and shift a list of images
    Rotate input output rotation
    Angle of rotation of the image in degrees. Positive angles will rotate the image counter-clockwise from the x axis.

    In the case of SAO 1-m longslit data, 180 degree rotation is required for easy wavelength calibration.
    """
    import glob
    import os, sys
    from pyraf import iraf
    iraf.images()
    iraf.imgeom()
    imlist = glob.glob(imlist_name)
    imlist.sort()   
    for i in range(len(imlist)):
        inim = imlist[i]
        print('Rotate '+inim+' with angle of '+str(round(rotation,3)) )
        iraf.rotate(input=inim, output='r'+inim, rotation=rotation)
        print('r'+inim+' is created.')

def sflip(imlist_name):
    """
    아크 이미지랑 사이언스 이미지랑 파장 방향이 다를때, 스펙트럼을 뒤집는다.
    세로로 되어 있는 이미지는 뒤집히나마나이므로, 가로방향으로 나오는 1d spectrum을 뒤집으면 된다.

    sflip -- Flip data and/or dispersion coordinates in spectra
    sflip input output

    """
    import glob
    import os, sys
    from pyraf import iraf
    iraf.noao()
    iraf.onedspec()
    imlist = glob.glob(imlist_name)
    imlist.sort()   
    for i in range(len(imlist)):
        inim = imlist[i]
        print('Flip spectrum')
        iraf.sflip(input=inim, output='r'+inim, coord_flip='yes', data_flip='yes')

def identify(imlist_name):
    """
    identify arc image for wavelength calibration.
    image section "column" is good.
    fitting function : "chebyshev", "legendre", "spline1", or "spline3"
    """
    import glob
    import os, sys
    from pyraf import iraf
    iraf.noao()
    iraf.imred()
    iraf.kpnoslit()
    imlist = glob.glob(imlist_name)
    imlist.sort()   
    for i in range(len(imlist)):
        inim = imlist[i]   
        print('Identification ongoing...')
        iraf.identify(images=inim, section='middle column', function='chebyshev', order=15, niterate=20, low_reject=3, high_reject=3 )

def hedit(imlist_name, wavfile):
    """
    Put the result of wavelenght calibration on standard, science images.
    hedit fcdbstd.ms.fits REFSPEC1 "fcdbwav.ms.fits" add+ ver- show+
    """
    import glob
    import os, sys
    from pyraf import iraf    
    imlist = glob.glob(imlist_name)
    imlist.sort() 
    for i in range(len(imlist)):
        inim = imlist[i]     
        iraf.hedit(images=inim, fields='REFSPEC1', value=wavfile, add='yes', verify='no', show='yes' )

def set_observatory(obsid):
    """
    obsid -- Examine and set observatory parameters

    Custom obsdb.dat is at /home/lim9/anaconda3/lib/python3.7/site-packages/lgpy/SAO_KL400

    """
    import glob
    import os, sys
    from pyraf import iraf
    iraf.noao()
    iraf.reset(obsdb='/home/lim9/anaconda3/lib/python3.7/site-packages/lgpy/SAO_KL400/obsdb.dat')
    iraf.observatory(command="set", obsid=obsid)


def setairmass(inim, obs, ra, dec):
    """
    setairmass -- update image headers with the effective airmass
    obs   -- SAO1-m : "sao" 
    ra    -- target ra hh:mm:ss
    dec   -- target dec dd:mm:ss
    location -- Put your observatory information.
    *** Hourly airmass for Feige 56 ***

    Epoch 2000.00: RA  12 06 47.3, dec  11 40 13
    Epoch 2019.35: RA  12 07 46.6, dec  11 33 45

    At midnight: UT date 2019 May 10, Moon 0.33 illum,  59 degr from obj

    Local      UT      LMST      HA     secz   par.angl. SunAlt MoonAlt

    22 00    11 00    10 40    -1 28    1.186   -33.6     -6.3    52.1
    23 00    12 00    11 40    -0 28    1.119   -12.4    -16.3    40.7
    0 00    13 00    12 40     0 32    1.121    14.2     ...     29.1
    1 00    14 00    13 40     1 32    1.194    34.7     ...     17.7
    2 00    15 00    14 40     2 33    1.363    46.2     ...      6.6
    3 00    16 00    15 40     3 33    1.701    51.8     ...     ... 
    4 00    17 00    16 41     4 33    2.436    53.9     ...     ... 
    5 00    18 00    17 41     5 33    4.691    53.8     ...     ... 
    6 00    19 00    18 41     6 33  (v.low)    51.8    -15.8    ... 
    7 00    20 00    19 41     7 33   (down)    47.9     -5.7    ... 
    """
    import glob
    import os, sys
    from pyraf import iraf
    from astropy.time import Time
    from astropy.io import fits
    from astropy import units as u
    hdr = fits.getheader(inim)
    t = Time(hdr['date-obs'], format='isot', scale='utc', location=(126.95333299999999*u.deg, 37.45704167*u.deg, 190*u.m))
    st = t.sidereal_time('apparent') # in hourangle
    set_observatory(obs)
    st_hms = str(int(st.hms.h)).zfill(2)+':'+str(int(st.hms.m)).zfill(2)+':'+str(int(st.hms.s)).zfill(2)
    iraf.hedit(images=inim, fields='st', value=st_hms, add='yes', verify='yes')
    iraf.hedit(images=inim, fields='ra', value=ra, add='yes', verify='yes')
    iraf.hedit(images=inim, fields='dec', value=dec, add='yes')
    print('OK.')
    iraf.astutil()
    iraf.setairmass(images=inim, observatory=obs, equinox='epoch', date = "date-obs", exposure = "exptime", airmass = "airmass",show = 'yes',ut="date-obs", override = 'yes')

def dispcor(imlist_name, database='database'):
    """
    dispcor -- Dispersion correct and resample spectra
    dispcor input output [records]
    database -- path in which idXXX file. ex) /data1/SN2019ein/work/SAO_Spectrum/red/20190509/arc/database/
     dispcor fcdbstd.ms.fits wfcdbstd.ms.fits
    """
    import glob
    import os, sys
    from pyraf import iraf   
    iraf.noao()
    iraf.imred()
    iraf.kpnoslit()
    imlist = glob.glob(imlist_name)
    imlist.sort() 
    for i in range(len(imlist)):
        inim = imlist[i]         
        iraf.dispcor(input=inim, output='w'+inim, database='database', linearize = 'no')
        iraf.splot(images='w'+inim)

def standard(imlist_name, obj, obs):
    """
    standard -- Add standard stars to sensitivity file
    standard input [records] output
    standard wfcdbstd.ms.fits (no) "sao"

    extinct = "/iraf/iraf/noao/lib/onedstds/iidscal/feige56.dat"

    caldir  = "/iraf/iraf/noao/lib/onedstds/iidscal/"
    """
    import glob
    import os, sys
    from pyraf import iraf   
    iraf.noao()
    iraf.onedspec()
    imlist = glob.glob(imlist_name)
    imlist.sort() 
    for i in range(len(imlist)):
        inim = imlist[i]  
        iraf.standard(input=inim, output='s'+inim[:-5], extinct='/iraf/iraf/noao/lib/onedstds/iidscal/feige56.dat', observatory=obs, caldir="/iraf/iraf/noao/lib/onedstds/iidscal/", star_name=obj) 

def sensfunc(standards, obs) :
    """
    sensfunc -- Determine sensitivity and extinction functions
    sensfunc standards sensitivity
    """
    import glob
    import os, sys
    from pyraf import iraf
    iraf.noao()
    iraf.onedspec()
    output_sens = "sens"
    iraf.sensfunc(standards=standards, sensitivity=output_sens, extinct='/iraf/iraf/noao/lib/onedstds/ctioextinct.dat', observatory=obs)    
    iraf.splot(output_sens)

def calibrate(imlist_name, obs) :
    """
    calibrate -- Apply extinction corrections and flux calibrations
    calibrate input output [records]
    """
    import glob
    import os, sys
    from pyraf import iraf   
    iraf.noao()
    iraf.onedspec()
    imlist = glob.glob(imlist_name)
    imlist.sort() 
    for i in range(len(imlist)):
        inim = imlist[i]  
        #iraf.calibrate(input=inim , output='s'+inim, extinct='/iraf/iraf/noao/lib/onedstds/ctioextinct.dat', flux='yes', observatory=obs)
        iraf.calibrate(input=inim , output='s'+inim, extinct='yes', flux='yes', observatory=obs)

def listpix(inim):
    """
    Extract calibrated (wavelength, flux) 1d spectrum into data file.
    Header information of object, exptime, date-obs(UT) are required to change file name terminally. 
    """
    import glob
    import os, sys
    from pyraf import iraf
    from astropy.io import fits
    from astropy.time import Time

    hdr     = fits.getheader(inim)
    obj     = hdr['object']
    exptime = hdr['exptime']
    obsdate = hdr['date-obs']
    t       = Time(obsdate)
    ccd     = 'SAO_KL400'
    output_data = 'Calib-'+ccd+'-'+obj+'-'+str(t.datetime.year).zfill(2)+str(t.datetime.month).zfill(2)+str(t.datetime.day).zfill(2)+'-'+str(t.datetime.hour).zfill(2)+str(t.datetime.minute).zfill(2)+str(t.datetime.second).zfill(2)+'-'+str(int(exptime))+'.data'
    iraf.imutil()
    org_stdout = sys.stdout
    f = open(output_data, "w+")
    sys.stdout = f
    iraf.listpixels(images=inim)
    sys.stdout = org_stdout
    f.close()
    os.system('pluma '+output_data+' &')

def fig(input_data) :
    """
    """
    import glob
    import os, sys
    from astropy.io import fits
    from astropy.io import ascii
    import matplotlib.pyplot as plt
    data = ascii.read(input_data)
    w = data['col1']
    f = data['col2']
    plt.plot(w, f)
    plt.show()
    plt.close()

def extdata(inim, plot=True): 
    """
    """
    import glob
    import os, sys
    import numpy as np
    from astropy.io import fits
    from astropy.io import ascii
    import matplotlib.pyplot as plt    
    data, hdr = fits.getdata(inim, header=True)
    flux = np.array(data)
    wave = np.ones(hdr['NAXIS1'], dtype=float)
    for i in range(hdr['NAXIS1']) :
        wave[i] = hdr['CRVAL1'] + i*hdr['CDELT1']




    
    
    
