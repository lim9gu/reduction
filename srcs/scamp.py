def projedit(imlist_name) :
    from pyraf import iraf
    from astropy.io import fits
    import glob
    imlist      = glob.glob(imlist_name)
    imlist.sort()    
    for i in range(len(imlist)):
        inim = imlist[i]
        data, hdr = fits.getdata(inim, header=True)
        try :
            del hdr['PV2_1']
        except:
            pass
        try :    
            del hdr['PV2_2']
        except:
            pass
        try :
            del hdr['PV2_3']
        except:
            pass
        #hdr['PROJP1'] = 1.0
        #hdr['PROJP3'] = -50
        del hdr['PROJP1']
        del hdr['PROJP3']
        del hdr['PROJP5']
        hdr['EQUINOX'] = '2000.0'
        hdr['RADECSYS'] = 'FK5'
        fits.writeto('p'+inim, data, header=hdr, overwrite=True)

def run_scamp(imlist_name, ccd, refcat='USNO-B1') :
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
    from lgpy.hdrcheck import wcscenter

    # UKIRT TAN-SIP
    if ccd == 'UKIRT_TAN' :
        cd11   =      -1.0251200248E-06 #/ Transformation matrix
        cd12   =       0.000111444818303 #/ no comment
        cd21   =      -0.000112119992813 #/ no comment 
        cd22   =      -3.79779999189E-07 #/ no comment 
    # UKIRT RA--ZPN DEC--ZPN
    elif ccd == 'UKIRT_ZPN' :
        cd11   =      -5.5747965E-05 #/ Transformation matrix
        cd12   =       9.9499104E-08 #/ no comment
        cd21   =      -8.9013383E-08 #/ no comment 
        cd22   =      -5.5714834E-05 #/ no comment  
    # SQUEAN
    elif ccd == 'SQUEAN' :
        cd11   = -7.9607532911290E-05 #/ Linear projection matrix              
        cd12   =  -1.754777100142E-06 #/ Linear projection matrix               
        cd21   =   -1.89497183878E-06 #/ Linear projection matrix               
        cd22   =   7.969879099597E-05 #/ Linear projection matrix 

    elif ccd == 'KCT_STX16803' : 
        cd11   = -0.000200202543613 #/ Linear projection matrix              
        cd12   =  1.91058495346E-06 #/ Linear projection matrix               
        cd21   =   -8.08992122764E-0 #/ Linear projection matrix               
        cd22   =   1.91058495346E-06 #/ Linear projection matrix     
    sexconfig   = '/data1/code/astrom/astrom.sex'
    sexconv     = '/data1/code/astrom/astrom.conv'
    sexnnw      = '/data1/code/astrom/astrom.nnw'
    sexparam    = '/data1/code/astrom/astrom.param'
    psfconfig   = '/data1/code/psfex.config/prepsfex.sex'
    scampconfig = '/data1/code/astrom/astrom.scamp'

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
        hdr = fits.getheader(inim)
        #ra, dec    = wcscenter(inim)
        ra = hdr['ra']
        dec = hdr['dec']
        #ra= '208.3658234'
        #dec= '40.25456995'
        #ra = 2.081642579100E+02
        #dec= 4.051871459365E+01
        rad   = coord.Angle(ra,unit=u.deg)    
        radd  = rad.degree
        decd  = coord.Angle(dec,unit=u.deg)
        decdd = decd.degree    
        xpix  = hdr['NAXIS1']/2.
        ypix  = hdr['NAXIS2']/2.
        #xpix  = 364.16726
        #ypix  = 361.08691
        #xpix  = 6.0251172E+03
        #pix  = -1.9102533E+03
        threshold,minarea = 3,3
        output_cat = 'astromtest.cat'
        sexcom = 'sex -c '+sexconfig +' '+inim+' -PARAMETERS_NAME '+sexparam+' -DETECT_THRESH '+str(threshold) +' -DETECT_MINAREA '+str(minarea)+' -CATALOG_NAME '+output_cat+' -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+' -CATALOG_TYPE FITS_LDAC'

        scampcom='scamp -c '+scampconfig+' astromtest.cat '+ '-ASTREF_CATALOG '+refcat+' -CHECKPLOT_DEV NULL -CROSSID_RADIUS 3.0 -PROJECTION_TYPE TAN'
        print(sexcom)
        print(scampcom)

        # iraf setting
        wcsreset(inim,'physical')
        wcsreset(inim,'world')
        hedit(inim,'WAT0_001','system=image')
        hedit(inim,'WAT1_001', 'wtype=tan axtype=ra')
        hedit(inim,'WAT2_001', 'wtype=tan axtype=dec')
        hedit(inim,'RADECSYS', 'FK5')
        hedit(inim,'EQUINOX', '2000.0')
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