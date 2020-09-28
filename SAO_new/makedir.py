def makedir(camera='STX16803'):
    """
    1. Description
    : Make directory of camera name (STX16803 or Kepler)

    2. Usage
    >>> 

    3. History
    2019.04.21 : G. Lim created.
    """
    import glob
    import os, sys
    workdir = '/data3/IMSNG/IMSNGgalaxies'
    curdir = os.getcwd()
    if curdir != workdir :
        print('Please run this code in '+workdir)
    elif curdir == workdir :
        folder = glob.glob('*')
        folder.sort()
        for i in range(len(folder)):
            os.chdir(folder[i])
            obs = glob.glob('*')
            obs.sort()
            if 'SAO' in obs :
                print('No camera dir. Make '+camera)
                os.chdir('SAO')
                os.system('/usr/bin/mkdir '+'./'+ camera)
                os.system('/usr/bin/mv Cal*.fits '+camera)
                #os.system('rsync -a B '+camera+'/B')
                os.chdir('../../')
            else :
                print(folder[i]+' has no SAO data.')
                os.chdir('../')
                pass
        print('Done.')