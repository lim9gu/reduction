def imcopy(imlist_name, region) :
    """
    SAO KL400 SN2019ein 
    [893:1083,337:1695]
    """
    import glob
    from pyraf import iraf
    imlist = glob.glob(imlist_name)
    imlist.sort()
    for i in range(len(imlist)):
        inim = imlist[i]
        iraf.imcopy(input=inim+region,output='t'+inim,verbose = 'yes')

def imcombine(imlist_name, newimage, combine= 'median', reject='sigclip', scale='none', zero='mode'):
    from pyraf import iraf
    import os, sys
    import glob
    image_to_com = glob.glob(imlist_name)
    iraf.imcombine(','.join(image_to_com), newimage, combine = combine, reject=reject, scale=scale, zero=zero)


