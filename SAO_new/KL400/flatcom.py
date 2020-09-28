def imarith(operand1, op, operand2, result):
    from pyraf import iraf
    print(operand1+' '+op+' '+operand2+' = '+result)
    iraf.imarith(operand1=operand1, op=op, operand2=operand2, result=result)

def flatcom(inim_list, outim_list, combine='median',reject='sigclip', scale='mode') :
    import os, sys
    from pyraf import iraf
    iraf.noao()
    iraf.imred()
    iraf.ccdred()
    iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
    iraf.flatcombine(input=inim_list, output=outim_list, combine=combine, reject=reject, process='no', scale=scale, ccdtype='' )
    print('Output masterflat is created.')    

def runflatcom(imlist_name, exptime) :
    import glob
    import sys, os
    import numpy as np
    from pyraf import iraf
    from astropy.io import fits
    import lgpy.hselect as hs
    imlist = hs.hselect(imlist_name, 'IMAGETYP', 'Flat Field')
    #imlist  = glob.glob(imlist_name)
    imlist.sort()
    print('Process ongoing...')
    input_list  = ','.join(imlist)
    output_list = 'flat.fits'
    flatcom(input_list, output_list, combine='median', reject='sigclip', scale='mode' )
    imarith(output_list, '-', 'zero.fits', 'z'+output_list )
    imarith('z'+output_list, '-', 'zdark'+str(int(exptime))+'.fits', 'd'+'z'+output_list )
    data, newhdr = fits.getdata('dz'+output_list, header=True)
    x = np.mean(data)
    nimage = data/x
    newflat_name = 'n'+output_list
    fits.writeto(newflat_name, nimage, header=newhdr, overwrite=True)

