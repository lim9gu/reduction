def biascom(inim_list, outim_list, combine='median',reject='minmax', scale='none'):
    import os, sys
    from pyraf import iraf
    iraf.noao()
    iraf.imred()
    iraf.ccdred()
    iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
    iraf.zerocombine(input=inim_list, output=outim_list, combine=combine, reject=reject, process='no', scale=scale, ccdtype='' )
    print('Output masterbias is created.')

def runbiascom(imlist_name) :
    import glob
    import sys, os
    from pyraf import iraf
    from astropy.io import fits
    from lgpy.hselect import hselect as hs

    imlist = hs(imlist_name, 'IMAGETYP', 'Bias Frame')
    #imlist = glob.glob(imlist_name)
    imlist.sort()
    #input_list = ','.join(imlist)
    #output_list = 'z'+',z'.join(imlist)
    print('Process ongoing...')
    def split_list(a_list):
        half = int(len(a_list)/2)
        return a_list[:half], a_list[half:]
    if len(imlist)%2 == 1 :
        input_list = ','.join(imlist)
        output_list = 'zero.fits'
        print('Odd number images, stack using median.')
        biascom(input_list, output_list, combine='median', reject='sigclip', scale='none')    
    elif len(imlist)%2 == 0 :
        inim_list1, inim_list2 = split_list(imlist)
        input_list1 = ','.join(inim_list1)
        input_list2 = ','.join(inim_list2)
        output_list1 = 'zero1.fits'
        output_list2 = 'zero2.fits'
        print('Even number images, stack twice (median) --> (average)')
        biascom(input_list1, output_list1, combine='median', reject='sigclip', scale='none')
        biascom(input_list2, output_list2, combine='median', reject='sigclip', scale='none')
        input_list = [output_list1, output_list2]

        output_list = 'zero.fits'
        biascom(','.join(input_list), output_list, combine='average', reject='sigclip', scale='none')


