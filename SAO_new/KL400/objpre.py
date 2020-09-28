def imarith(operand1, op, operand2, result):
    from pyraf import iraf
    print(operand1+' '+op+' '+operand2+' = '+result)
    iraf.imarith(operand1=operand1, op=op, operand2=operand2, result=result)

def objpre(imlist_name):
    """
    bias subtraction
    dark subtraction
    flat fielding
    """
    import glob
    import os, sys
    from astropy.io import fits
    imlist  = glob.glob(imlist_name)
    imlist.sort()
    print('Process ongoing...')    
    for i in range(len(imlist)) :
        inim = imlist[i]
        data, hdr = fits.getdata(inim, header=True)
        exptime = hdr['exptime']
        print('bias subtraction...')
        imarith(inim, '-', 'zero.fits', 'z'+inim)
        print('dark subtraction...')
        imarith('z'+inim, '-', 'zdark'+str(int(exptime))+'.fits', 'dz'+inim)
        print('flat fielding...')
        imarith('dz'+inim,'/', 'nflat.fits', 'fdz'+inim)