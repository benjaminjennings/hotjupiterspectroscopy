#---------------------------------------------------------------------------------------------

#Benjamin Jennings
#16317714

#---------------------------------------------------------------------------------------------

#------------SCRIPT TO PRODUCE MODEL TRANSMISSION SPECTRA TO ANALYSE ATMOSPHERE---------------

#---------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import h5py

#---------------------------------------------------------------------------------------------

#MAKE SURE TO BE IN THE RIGHT DIRECTORY FOR THE TEMPLATE YOU ARE CREATING

#---------------------------------------------------------------------------------------------

#function to List the Species available to create templates from for a given temperature tmp

def LS(tmp):

    """Returns list of available species for a given Temperature string (eg: '2000K')."""

    with h5py.File('/cphys/ugrad/2017-18/SF/BEJENNIN/Desktop/Project/My_Templates/CrossSections.h5','r') as F:
        
        #-------------------------------------------------------------------------------------

        #print available cross-sections
        
        print('')
        print('Temperatures: {}'.format(F.keys()))
        print('')
        print(F['{}'.format(tmp)].keys())
        print(len(F['{}'.format(tmp)].keys()))
        
        #-------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------

#defining my Template Generator function 

def TG(tmp,sp):
    
    """Returns a template for a given species and temperature"""
    
    with h5py.File('/cphys/project/BEJENNIN_project/Project/My_Templates/CrossSections.h5','r') as F:
        
        #-------------------------------------------------------------------------------------

        #print available cross-sections
        #print F.keys()
        #print F['{}'.format(tmp)].keys() 
    
        #refractive index is wavelength dependent (and also depends on the medium slightly)
        #import all for key
        w = F['w_air'][:] 
    
    
        model = F['{}/spec_{}'.format(tmp,sp)][:]  / 1.e4 
        
        #-------------------------------------------------------------------------------------
        
        #constants in SI units
        
        RSun = 695000000. #m
        RJup = 69911000.
        MJup = 1.9e+27
        G = 6.6742e-11
        mh = 1.67e-27 #kg
        k = 1.3806503e-23
        
        Rstar = 1.458*RSun 
        Mpl = 1.184*MJup
        Rpl = 0.118*Rstar
        g = G * Mpl / Rpl**2
        bar = 100000. #Pa
        T = 1500 #K
        
        #scale height
        Hs = k * T / (2.3*mh * g)
        
        #compute the constant C.
        gamma = 0.56
        P0 = 1.e1 * 100000. #reference pressure in bar - ie pressure at Rpl
        C = Rpl + Hs * (gamma + np.log(P0 / g) + 0.5*np.log(2.*np.pi*Rpl/Hs) - np.log(2.3*mh))
        
        #finally get the transmission spectrum
        mix = 1e-6 #mixing ratio of species
        R = C + Hs * np.log(mix * model)

        #-------------------------------------------------------------------------------------
        
        #set pressure of cloud deck
        P_cloud = 1e-3 * 100000. #in bar

        R_cloud = Hs * np.log(P0/P_cloud) + Rpl

        #truncate the continuum model with cloud deck
        R_cont = np.maximum(R,R_cloud)
        
        #-------------------------------------------------------------------------------------

        #planet to star radius ratio
        RpRs = R_cont / Rstar

        #transit depth
        F = RpRs ** 2

        #or with subtracted continuum of cloud
        dF = RpRs ** 2 - (R_cloud/Rstar) ** 2
        
        #-------------------------------------------------------------------------------------

        #save the transmission spectrum (in negative flux units)

        #filename = 'W121b_spec_{}.npy'.format(sp)
        filename = 'W121b_spec_{}_{}.npy'.format(sp,tmp)

        np.save(filename,np.array([w,-dF]))
        
        #-------------------------------------------------------------------------------------
        
        #plotting the resulting template
        
        #plt.figure(figsize=(12,6))
        
        #plt.title('{} Template ({})'.format(sp,tmp))
        #plt.xlabel('$\lambda$ ($cm$)')
        #plt.plot(w,-dF,'maroon')
        #plt.show()
        
        #-------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------

#template plotting function

def Tplot(tmp,sp):
    
    #loading the required template
    
    b = np.load('/cphys/ugrad/2017-18/SF/BEJENNIN/Desktop/Project/My_Templates/Generated_Templates/W121b_spec_{}.npy'.format(sp))

    w_model,model = b[0],b[1]
    
    #plotting the resulting template
        
    fig = plt.figure(figsize=(12,6))
        
    plt.title('{} Template ({})'.format(sp.replace('_',''),tmp))
    plt.xlabel('$\lambda$ ($\AA$)')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.plot(w_model,model,'maroon')
    plt.show()

#---------------------------------------------------------------------------------------------
