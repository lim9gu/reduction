#/home/lim9/anaconda3/lib/python3.7/site-packages/lgpy/SAO_KL400/darkcom.py
def imarith(operand1, op, operand2, result):
    from pyraf import iraf
    print(operand1+' '+op+' '+operand2+' = '+result)
    iraf.imarith(operand1=operand1, op=op, operand2=operand2, result=result)

def explist(imlist_name):
    import glob
    from astropy.io import fits
    if type(imlist_name) == str :
        dark = glob.glob(imlist_name)
    elif type(imlist_name) == list :
        dark = imlist_name
    allexptime = []
    for i in range(len(dark)) :
        hdr = fits.getheader(dark[i])
        allexptime.append(hdr['exptime'])
    expset = set(allexptime)
    exptime = list(sorted(expset))
    return exptime

def darkcom(inim_list, outim_list, combine='median',reject='minmax', scale='none'):
    import os, sys
    from pyraf import iraf
    iraf.noao()
    iraf.imred()
    iraf.ccdred()
    iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
    iraf.darkcombine(input=inim_list, output=outim_list, combine=combine, reject=reject, process='no', scale=scale, ccdtype='' )
    print('Output masterdark is created.')

def split_list(a_list):
    half = int(len(a_list)/2)
    return a_list[:half], a_list[half:]

def rundarkcom(imlist_name,  biassub = True) :
    import glob
    import sys, os
    import numpy as np
    from pyraf import iraf
    from astropy.io import fits
    import lgpy.hselect as hs
    dark = hs.hselect(imlist_name, 'IMAGETYP', 'Dark Frame')
    exptime = explist(dark)
    for i in range(len(exptime)) :
        exp = int(exptime[i])
        imlist = hs.hselect(imlist_name, 'exptime', exp)
        imlist.sort()
        print('Process ongoing...')
        if len(imlist)%2 == 1 :
            input_list = ','.join(imlist)
            output_list = 'dark'+str(exp)+'.fits'
            print('Odd number images, stack using median.')
            darkcom(input_list, output_list, combine='median', reject='sigclip', scale='none')   
        elif len(imlist)%2 == 0 :
            inim_list1, inim_list2 = split_list(imlist)
            input_list1 = ','.join(inim_list1)
            input_list2 = ','.join(inim_list2)
            output_list1 = 'dark'+str(exp)+'.1.fits'
            output_list2 = 'dark'+str(exp)+'.2.fits'
            print('Even number images, stack twice (median) --> (average)')
            darkcom(input_list1, output_list1, combine='median', reject='sigclip', scale='none')
            darkcom(input_list2, output_list2, combine='median', reject='sigclip', scale='none')
            input_list = [output_list1, output_list2]
            output_list = 'dark'+str(exp)+'.fits'
            darkcom(','.join(input_list), output_list, combine='average', reject='sigclip', scale='none')
        if biassub == True:
            try :
                imarith(output_list, '-', 'zero.fits', 'z'+output_list)
            except :
                print('No masterbias image or unusable file... ERROR.')


