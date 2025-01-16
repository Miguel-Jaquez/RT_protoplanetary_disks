### plot polarization results
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

def plot_test(image_stokes, image_tau, figure_to_plot='PI', output_model='', nx=10, ny=10):
    #plt.figure()
    fig, ax = plt.subplots(1,1, figsize=(10,10))
    im_tau = image.readImage(image_tau)
    data_tau = im_tau.image[:,:,0]

    im     = image.readImage(image_stokes)
    dataI  = im.image[:,:,0,0]
    dataQ  = im.image[:,:,1,0]
    dataU  = im.image[:,:,2,0]
    dataPI = np.sqrt(im.image[:, :, 1,0]**2 + im.image[:, :, 2,0]**2 ) # Polarized intensity
    dataP  = np.sqrt(im.image[:, :, 1,0]**2 + im.image[:, :, 2,0]**2 ) / im.image[:, :, 0,0]    #polarization fraction

    #
    if figure_to_plot == 'I':
        dataI[dataI<= 1e-20] = 0
        data = dataI/np.nanmax(dataI)
        #print ('min_data', data.mean())
        #cb_label = 'I' + r'$_\nu$' + ' [erg/s/cm/cm/Hz/ster]'
        cb_label = 'I' + r'$_\nu$/'+ 'I' + r'$_{\nu, max}$'
        #out_name = 'I_{}_edgeon_{}_disk'.format(model_output,type_disk)
    elif figure_to_plot == 'P':
        data = dataP * 100
        print (np.nanmax(data))
        cb_label = 'p' + r'$_f$' + '[%]'
        #vmin = 0
        #vmax = 1
        #out_name = 'Pf_{}_edgeon_{}_disk'.format(model_output,type_disk)
    elif figure_to_plot == 'PI':
        data = dataPI/dataPI.max()
        #cb_label = 'PI' + r'$_\nu$' + ' [erg/s/cm/cm/Hz/ster]'
        cb_label = 'PI' + r'$_\nu/$'+ 'PI' + r'$_{\nu, max}$'
        #out_name = 'PI_{}_edgeon_{}_disk'.format(model_output, type_disk)

    # Select the coordinates of the data
    x = im.x / natconst.au
    #print (im.x)
    y = im.y / natconst.au
    xlab = 'X [au]'
    ylab = 'Y [au]'

    ext = (x[0], x[im.nx - 1], y[0], y[im.ny - 1]) #sacar esto del header
    #print(x[im.nx - 1])

    # Now finally put everything together and plot the data
    #plot tau contours 
    c = ax.contour(data_tau.T, levels=5,
                   #levels = [1, 3,5,7, 10],  
                   colors='w', linestyles='dashed', extent=ext)
    ax.clabel(c, inline=1, fontsize=20)

    c_map = plt.cm.jet
    c_map.set_bad('w')
    data_plot = data.T.copy()
    data_plot[data_plot<1e-10] = np.nan

    if figure_to_plot == 'P': # lo dejamos con vmax=1 que ninguno pasa de 1 % de polarizacion
        implot = ax.imshow(data_plot, extent=ext, cmap='hot_r') #, norm=LogNorm(vmin=0.001, vmax=1))# vmin= 0, vmax=0.8) #, interpolation=interpolation, **kwargs)

    else:
        implot = ax.imshow(data_plot, extent=ext, cmap=c_map, norm=LogNorm(vmin=0.001, vmax=1)) #, interpolation=interpolation, **kwargs)

    '''
    #set the columns names
    if (i==0 or i == 1 or i==2) and 'inclination' in model_output: # top row put the names of the columns
        print (i)
        #ax.set_title(r'$\lambda$=' + ("%.5g" % im.wav[ifreq]) + r'$\mu$m')
        ax.set_title(inclinations[i])
    elif (i==0 or i == 1 or i==2) and 'grain_size' in model_output:
        print (i)
        ax.set_title(r'$a_{max}$=' + ("%.5g" % grain_sizes[i] + r'$\mu$m'))
    else:
        pass
    # set the row names
    if (i==0 or i ==3 or i==6 or i==9 or i==12 or i==15) and 'inclination' in model_output:
        #print (i)
        wavelength_text = ("%.5g" % im.wav[ifreq]) + r'$\mu$m'
        ax.text(-0.3,0.35, wavelength_text ,color='k', rotation='vertical',fontsize=15, transform=ax.transAxes)
    elif (i==0 or i ==3 or i==6 or i==9 or i==12 or i==15) and 'grain_size' in model_output:
        wavelength_text = ("%.5g" % im.wav[ifreq]) + r'$\mu$m'
        ax.text(-0.3,0.35, wavelength_text ,color='k', rotation='vertical',fontsize=15, transform=ax.transAxes)
    else:
        pass
    '''
    cbar = plt.colorbar(implot)
    cbar.set_label(cb_label)

    ########### plot the polarization direction


    # ext = (x[0], x[image.nx-1], y[0], y[image.ny-1])
    iix = [int(np.floor(i)) for i in np.arange(nx) * float(x.shape[0]) / nx]
    iiy = [int(np.floor(i)) for i in np.arange(ny) * float(x.shape[0]) / ny]
    xr = x[iix]
    yr = y[iiy]
    xxr, yyr = np.meshgrid(xr, yr, indexing='ij')
    grid_dum = np.meshgrid(xr, yr, indexing='ij')

    #mask_plot = dataI==0
    #dataQ[mask_plot] = np.nan
    #dataU[mask_plot] = np.nan
    qqr = (np.squeeze(dataQ)/np.squeeze(dataI).clip(1e-60))[np.ix_(iix, iiy)]
    uur = (np.squeeze(dataU)/np.squeeze(dataI).clip(1e-60))[np.ix_(iix, iiy)]
    #qqr = (np.squeeze(im.image[:, :, 1, ifreq])
     #      / np.squeeze(im.image[:, :, 0, ifreq]).clip(1e-60))[np.ix_(iix, iiy)]
    #uur = (np.squeeze(im.image[:, :, 2, ifreq])
     #      / np.squeeze(im.image[:, :, 0, ifreq]).clip(1e-60))[np.ix_(iix, iiy)]
    lpol = np.sqrt(qqr**2 + uur**2).clip(1e-60)
    qqr /= lpol
    uur /= lpol
    ang = np.arccos(qqr) / 2.0
    ii = (uur < 0)
    if True in ii:
        ang[ii] = np.pi - ang[ii]
    vx = np.cos(ang)
    vy = np.sin(ang)
    ii = (lpol < 1e-6)
    vx[ii] = 1e-6 #0.0001
    vy[ii] = 1e-6 #0.0001
    # tamaño de las felchas dependiente del grado de polarizacion
    #degreeLength = 1/(lpol.flatten()*2.2)/max(float(len(xr))/figsize[0]/6, float(len(yr))/figsize[1]/3)

    ax.quiver(xxr, yyr, vx, vy, color='k', pivot='mid', 
            # tamaño de las flechas
            #scale=1/degreeLength, 
            #units='inches', scale_units='inches', angles='xy',
            #headwidth=0, headlength=1, headaxislength=1, width=0.02)
            scale=2. * np.max([nx, ny]), 
            #angles = 'xy',
            headwidth=1e-10, headlength=1e-10, headaxislength=1e-10)

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(output_model)
    #fig.supxlabel(xlab)
    #fig.supylabel(ylab)
    #fig.suptitle(model_output)
    plt.tight_layout()    
    plt.savefig(output_model)
    return ('finish plot')
    
def plot_tau(image_tau, output_model='', nx=10, ny=10):
    #fig_params = [1000, 45] #wav, incl
    #image_tau = "image_wav{}_incl{}_tau.out".format(fig_params[0], fig_params[1])
    #output_name = "image_wav{}_incl{}_tau.png".format(fig_params[0], fig_params[1])
    nx=10 
    ny=10
     
    plt.figure()
    fig, ax = plt.subplots(1,1)
    im_tau = image.readImage(image_tau)
    data_tau = im_tau.image[:,:,0]

    print (data_tau.max, data_tau.min)

    ### Select the coordinates of the data
    x = im_tau.x / natconst.au
    y = im_tau.y / natconst.au
    xlab = 'X [au]'
    ylab = 'Y [au]'

    ext = (x[0], x[im_tau.nx - 1], y[0], y[im_tau.ny - 1]) #sacar esto del header

    ### Now finally put everything together and plot the data
    ### plot tau contours 
    plt.imshow(data_tau.T, extent=ext)
    ### This is an attempt to automate the choice of levels.
    #tau_values = data_tau.flatten()
    #tau_values = tau_values[tau_values>0]
    #levels = [np.percentile(tau_values,100)] 

    c = ax.contour(data_tau.T, [0.001, 0.01, 0.1, 1.0, 5.0],  colors=['w','w','w','orange','r', 'orange'], linestyles='dashed', extent=ext)
    ax.clabel(c, inline=1, fontsize=10)

        
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    #ax.set_title("Optical depth at {} microns and inclination {}".format(wavelength, inclination))
    ax.set_title(output_model)
    #fig.supxlabel(xlab)
    #fig.supylabel(ylab)
    #fig.suptitle(model_output)
    plt.tight_layout()    
    plt.savefig(output_model)
    #plt.show()
    return print ('finish plot')



