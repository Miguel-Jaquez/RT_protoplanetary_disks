import numpy as np
import matplotlib.pyplot as plt
import os
import sys
#import astropy.units as asu
import subprocess
#sys.path.append('/fs/posgrado30/other0/jesus/radmc3dPy/lib/python3.9/site-packages')
sys.path.append('/opt/RADMC3D/radmc3d-2.0/python/radmc3dPy/radmc3dPy')
from radmc3dPy import *
from astropy.io import fits
# add radmc not neccesary in multivac 
radmc_path = '/opt/RADMC3D/radmc3d-2.0/src' #'/home/jesus/bin'
radmc_rl = os.path.join(radmc_path, 'radmc3d')
print (radmc_rl)
#the model and inp files are already created so we only need create the images

#create the image.out
bands = [870.0, 1000.0, 1300.0, 1500.0, 1800.0, 2100.0, 2500.0, 3100.0, 4000.0, 5000.0, 7000.0]
inclinations = [0.0,45.0]
for inclination in inclinations:
	for band in bands:
	    #image.makeImage(npix=300., wav=band, incl=90., phi=0., sizeau=350., stokes=True,  setthreads=80)
	    subprocess.run(radmc_rl+' image npix 300 lambda {} incl {} phi 0. sizeau 350. stokes setthreads 50'.format(band, inclination), shell=True, executable='/bin/bash')
	    os.system('mv image.out image_wav{}_incl{}_stokes.out'.format(int(band), int(inclination)))
	    #image.makeImage(npix=300, wav=band, incl=90., phi=0., sizeau=350., setthreads=50, tracetau=True)
	    subprocess.run(radmc_rl+' image npix 300 lambda {} incl {} phi 0. sizeau 350. tracetau setthreads 50'.format(band, inclination), shell=True, executable='/bin/bash')
	    os.system('mv image.out image_wav{}_incl{}_tau.out'.format(int(band), int(inclination)))

print ('finish')

