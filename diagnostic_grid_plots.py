### disgnostic plots

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
os.getcwd()
#%matplotlib inline
#%matplotlib widget

# General Properties of the Model
data = analyze.readData(ddens=True)#, dtemp=True)
print (data.taux)
total_dust_mass = data.getDustMass() 
print ('The total dust mass in the model is {:.2e} Msun'.format(total_dust_mass/natconst.ms))
# this is the total optical depth BUT per axis, eg x=R

data_temp = analyze.readData(dtemp=True)
data_temp.readDustTemp()
print ("temperature grid properties")
print (data_temp.dusttemp.min(), data_temp.dusttemp.mean(), data_temp.dusttemp.max())

# plot the dust properties
opac = analyze.readOpac(ext='optool_2.2295-0.02031-2.08_a100um', scatmat=True)

total_optical_depth = data.getTau(wav=0.55, axis='x') #Opacity at 0.55um :  20417.708186833515
# print the surface density
#surf_dens = data.getSigmaDust()
#print ('The surface density is ', surf_dens, 'g/cm**2')

print ("plot the dust properties")

#create theh layout
layout = [["kplot",'f11plot'],['f21/f11plot','fsplot']]

fig, axd = plt.subplot_mosaic(layout, figsize=(10,10))
analyze.plotDustOpac(opac, var = 'kabs', label = 'kabs', xlabel= 'Wavelength [$\mu$m]', ylabel= '$k_{abs}$, $k_{sca}$, $k_{ext}$ [$cm^{2} \, g^{-1}$]' ,ax=axd['kplot'])
analyze.plotDustOpac(opac, var = 'ksca', label = 'ksca', xlabel= 'Wavelength [$\mu$m]', ylabel= '$k_{abs}$, $k_{sca}$, $k_{ext}$ [$cm^{2} \, g^{-1}$]', ax=axd['kplot'])
axd['kplot'].set_xscale('log')
axd['kplot'].set_yscale('log')
axd['kplot'].legend()
#plt.xlim(1,1e4)
axd['kplot'].set_ylim(1,200)
axd['kplot'].set_title('Dust Properties')

#axd['f11plot'].set_title["Scattering "]
analyze.plotScatmat(opac, var='z11', idust=0, iwav=None, wav=1000, xvar='ang', iang=None, ang=None, ax=axd['f11plot'], xlabel= 'Angle [degree]', ylabel= 'z11', title=None)

analyze.plotScatmat(opac, var='linpol', idust=0, iwav=None, wav=1000, xvar='ang', iang=None, ang=None, ax=axd['f21/f11plot'], xlabel= 'Angle [degree]', title=None)

plt.savefig('Dust_properties.png')

print (data.taux.shape)

plt.figure()
# plot density: que estoy graficando ? o como mide la optical depth
c = plt.contourf(data.grid.x/natconst.au, np.pi/2. - data.grid.y, np.log10(data.rhodust[:,:,0,0].T), 10) #, vmin=-25, vmax=-16)
#colorbar(contourf_,ticks=range(vmin, vmax+3, 3))
cb = plt.colorbar(c ) #, ticks=range(-20, -14, 3) )#, norm = Norm())
cb.set_label(r'$\log_{10}{\rho} \, [\, g \, cm^{-3}]$')#, rotation=270.)
ct = plt.contour(data.grid.x/natconst.au, np.pi/2.-data.grid.y, data.taux[:,:,0].T, levels = [0.1,1.0,10], 
                 colors=['w','gray','k'], linestyles='solid') #,vmin=1e-10)
plt.clabel(ct, inline=1, fontsize=10)
plt.ylim(-0.2,0.2)
plt.xlabel('r [AU]')
plt.ylabel(r'$\pi/2-\theta$')
#plt.xscale('log')
plt.tight_layout()
plt.savefig("Density_and_opticaldepth_contours.png")


plt.figure()
c = plt.contourf(data_temp.grid.x/natconst.au, np.pi/2.-data_temp.grid.y, data_temp.dusttemp[:,:,0,0].T, [5,10,20,30,40,50,100,200])
#c = plt.contourf(data_temp.grid.x/natconst.au, np.pi/2.-data_temp.grid.y, data.rhodust[:,:,0,0].T, 30)

plt.xlabel('r [AU]')
plt.ylabel(r'$\pi/2-\theta$')
cb = plt.colorbar(c)
cb.set_label('T [K]', rotation=270.)

ct = plt.contour(data_temp.grid.x/natconst.au, np.pi/2.-data_temp.grid.y, data_temp.dusttemp[:,:,0,0].T, levels=[25,30,50,100,200],  colors='k', linestyles='solid')
plt.clabel(ct, inline=1, fontsize=10)
#plt.xscale('log')
plt.ylim(-0.5,0.5)
plt.savefig("Temperature_distribution_contours.png")





