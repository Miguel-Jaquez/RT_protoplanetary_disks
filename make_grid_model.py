#####-----------------------
# import the modules
#####-----------------------
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import astropy.units as asu
from matplotlib.colors import LogNorm
import os 
myhost = os.uname()[1]
if myhost == 'posgrado30.crya.privado':
    main_path = '/fs/posgrado30/other0/jesus/paper3/single_gs/'
    sys.path.append('/fs/posgrado30/other0/jesus/radmc3dPy/lib/python3.9/site-packages')

elif myhost == 'jesus-Latitude-3310':
    main_path = '/home/jesus/paper3/single_gs/'
    sys.path.append('/home/jesus/radmc3d-2.0/python/radmc3dPy/radmc3dPy')
from radmc3dPy import *
from astropy.io import fits
#os.getcwd()
#%matplotlib inline
#%matplotlib widget

print ('The  model is:')
print (models.getModelDesc('ppdisk'))

print ('Create the original parameters file')
## check if the problem_params.inp exist
if os.path.exists('./problem_params.inp'):
    print ("the file problem_params exist")
    #par = analyze.readParams()
    #par.printPar()
else: # if the file don exist is created
    analyze.writeDefaultParfile('ppdisk')

    ####-------
    # Modify the parameters file
    ####-------
    with open('problem_params.inp','r') as f:
        filedata = f.read()
        #f.close()
    #use the  exponential outer tapering model    
    filedata = filedata.replace('sigma_type                = 0','sigma_type                = 1') 
    # Rc    
    filedata = filedata.replace('hrpivot                   = 100.0*au','hrpivot                   = 50.0*au') 
    # Power exponent of the surface density distribution as a function of radius
    filedata = filedata.replace('plsig1                    = -1.0','plsig1                    = -0.2') 
    # Flaring index
    filedata = filedata.replace('plh                       = 1./7.','plh                       = 1.25') 
    # Ratio of the pressure scale height over radius at hrpivot
    filedata = filedata.replace('hrdisk                    = 0.1','hrdisk                    = 0.01') 
    # Outer radius of the disk
    filedata = filedata.replace('rdisk                     = 100.0*au','rdisk                     = 150.0*au') 
    # Inner radius of the disk
    filedata = filedata.replace('rin                       = 1.0*au','rin                       = 1.0*au') 
  
    with open('problem_params.inp','w') as f:
        f.write(filedata)
    f.close()
    #par = analyze.readParams()
    #par.printPar()
    print (filedata)

#####-------
# generate the gas distribution and create the .inp files from radmc
#####-------
setup.problemSetupDust('ppdisk', binary=False) #_with_temp_by_jaquez',writeDustTemp=True) #, writeGasTemp=True, binary=False)  #change binary to True to be faster 

####### ---------------
# Calculate the temperature as function of r for now, # ad function of phi will be added later as a general case 
#######
#read the grid properties
data = analyze.readData(ddens=True)#, dtemp=True)

#params
plht = -0.5
hrpivot = 50 * natconst.au   #ppar
R0 = 50 * natconst.au
#print ('{:.2}'.format(R0))
T0 = 30         #ppar
hrdisk = 0.01   #ppar

# the pass to cylyndical coordinates
rr, th = np.meshgrid(data.grid.x, data.grid.y)
zz = rr * np.cos(th)
rcyl = rr * np.sin(th)

# calculate the temperature as a function of r, phi
ht = np.zeros([data.grid.nx, data.grid.ny, data.grid.nz], dtype = np.float64)
#dum  = ppar['temp0'] * (rcyl/ppar['trpivot'])**ppar['plt'] # *rcyl?
dum  = T0 * (rcyl/R0)**plht # * rcyl

dum = dum.swapaxes(0,1)
for iz in range(data.grid.nz):
    ht[:,:,iz] = dum

# write the temperature file
#data = analyze.readData(ddens=True)#, dtemp=True)
fname = 'dust_temperature.dat'
wfile = open(fname, 'w')
hdr = np.array([1, data.grid.nx*data.grid.ny*data.grid.nz,1], dtype=int)
hdr.tofile(wfile, sep=" ", format="%d\n")
# Now we need to flatten the dust density array since the Ndarray.tofile function writes the 
# array always in C-order while we need Fortran-order to be written
data_temp_ = np.swapaxes(ht,0,2)
data_temp_.tofile(wfile, sep=" ", format="%.9e\n")
print ('finish write temperatura')
