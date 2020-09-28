def astroscrappy(imlist_name) :
    """
    1. Description 
    : Remove cosmic ray using astroscrappy which is fast proceeing LA cosmic. 

    2. Usage 
    >>> astroscrappy()

    3. History
    < 2016        C. Choi created.
    2019.01.26    G. Lim brought to unified processing pipeline. 
    2019.08.19    Edited for SAO KL400 version by G.Lim 
    """
    import os
    import glob
    import numpy as np
    import astroscrappy as cr
    import lgpy.obs as obs
    from astropy.io import fits

    os.nice(20)
    imlist = glob.glob(imlist_name)
    imlist.sort()
    for i in range(len(imlist)) :
        inim = imlist[i]
        def ascrapp(im, gain=obs.SAO_KL400.gain(), rdnoise=obs.SAO_KL400.rdnoise()):
            print im
            data, hdr = fits.getdata(im, header=True)
            c1,c2 = cr.detect_cosmics(data,gain=gain,readnoise=rdnoise, sigclip=4.5, sigfrac=0.3, objlim=5,niter=4,verbose=True)
            hdr.set('COMMENT','AstroSCRAPPY cr removed')
            fits.writeto('c'+im,c2,header=hdr)
        for im in lists : 
            ascrapp(im)
            os.system('ls c*.fits | wc -l && ls tbg*.fits |wc -l') 
    print 'astroscrappy is finished.'