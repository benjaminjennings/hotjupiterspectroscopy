#---------------------------------------------------------------------------------------------

#Benjamin Jennings
#16317714

#---------------------------------------------------------------------------------------------

#---------------FUNCTIONS FOR ANALYSING THE ATMOSPHERE OF HOT JUPITER EXOPLANETS--------------

#---------------------------------------------------------------------------------------------

import numpy as np 
import matplotlib.pyplot as plt

#ggplot parameters

plt.style.use('ggplot')
plt.rcParams['axes.facecolor']='w'
plt.rcParams['axes.edgecolor']='#999999'

#from numba import jit

#---------------------------------------------------------------------------------------------


#Function to perform the Cross-Correlation between a data set and a spectral template
#@jit(nopython=True)
def crosscorr(data,template):
    
    """
    Linear Cross-Correlation of Template and Data
    Need both parameters to be normalised to 0.
    """
    
    crc = np.empty((data.shape[0],data.shape[1]))

    for i in range(data.shape[0]):
        
        crc[i] = np.correlate(data[i],template,'same')

    return crc


#---------------------------------------------------------------------------------------------


#Function to test a range of velocities for the exoplanet.
#Shifts to planets rest frame and sums for each velocity.

def shiftsum(Kp,ph,vsys,crc):
    
    crcmap = np.empty((Kp.size,vsys.size))
    
    for j,k in enumerate(Kp):
        
        crcshift = np.empty(crc.shape)
        
        vp = k*np.sin(ph*2.*np.pi) 
        
        for i,dv in enumerate(vp): 
            
            crcshift[i] = np.interp(vsys+dv,vsys,crc[i])
  
        crcmap[j] = crcshift.sum(axis=0)
    
    return crcmap
    

#---------------------------------------------------------------------------------------------

    
#plotting function for the '3D' colormaps

def splot(x,y,data,c_map='bone',xlim1=-200,xlim2=200,ylim1=-500,ylim2=500,grid=0,cbar=0,x_title='Detection Significance',xlbl='$\mathrm{V_{sys}}$ $\mathrm{(km/s)}$',ylbl='$\mathrm{K_{p}}$ $\mathrm{(km/s)}$',show=1,save=0,closeplot=0,fn=''):
    
    """This function plots the shifted cross correlation for a range of velocities"""
    
    fig = plt.figure(figsize=(12,6))
    plt.pcolormesh(x,y,data,cmap=c_map)
    plt.title(x_title,fontsize='14',color='k')
    plt.xlabel(xlbl,fontsize='10',color='k')
    plt.ylabel(ylbl,fontsize='10',color='k')
    plt.xlim(xlim1,xlim2)
    plt.ylim(ylim1,ylim2)
    
    if grid == 1:
        
        mv = np.where(data == data.max())
        py1 = [y[int(mv[0])],y[int(mv[0])]]
        px1 = [xlim1,xlim2]
        
        px2 = [x[int(mv[1])],x[int(mv[1])]]
        py2 = [ylim1,ylim2]
        
        plt.plot(px1,py1,'k--',linewidth=0.5)
        plt.plot(px2,py2,'k--',linewidth=0.5)
        plt.text(xlim2-65,ylim2-10, '$V_{{sys}}$ = {}$km/s$, $K_{{p}}$ = {}$km/s$'.format(round(x[int(mv[1])],2),round(y[int(mv[0])],2)),bbox=dict(facecolor='grey',edgecolor='black',alpha=0.6))
    
    if cbar == 1:
        
        cbr = plt.colorbar(pad=0.01)
        cbr.set_label('$\mathrm{\sigma}$',rotation=0,fontsize=15,color='#838383')
       
    if save == 1:
    
        plt.savefig(fn,dpi=fig.dpi)
        
    if closeplot == 1:
        
        plt.close()
        
    else:
        
        plt.show()
    
#---------------------------------------------------------------------------------------------
    
    
    
    
    
