# -*- coding: utf-8 -*-
"""
Created on Fri May 19 17:38:34 2023

@author: alexn
"""

import numpy as np

#packages

def myQuadrupoleSignal(freq,params):
    """ Compute myQuadrupoleSignal time-domain waveform for compact binary coalescences.
        --------
        params dictionary with the following
        key-arguments :

        chipmass1 = chirpmass [Solarmass]
        d1 = luminosity distance [Mpc]
        iota = inclination angle [rad]
        phiRef = reference phase [rad]
        t_gps =  GPS trigger time [s]
        seglen = segment duration [s]
        --------
        t = array
        hp = array
        plus polarization
        hc = array
        cross polarization
    """
    G = 6.67430*(10**(-11))
    c = 299792458 
    #natural constants when working in SI-units
    #           
    chirpmass1   =   params['mchirp']       #luminosity chirpmass
    jota         =   params['iota']
    d1           =   params['distance']     #luminosity distance
    srate        =   params['srate']
    phiRef       =   params['phi_ref']
    t_gps        =   params['t_gps']
    seglen       =   params['seglen']
    time_shift   =   params['time_shift']
    #parameters
    
    #Definitions
    Solarmass = 1.98847*(10**(30))
    chirpmass = Solarmass*chirpmass1
    d = d1*3.0857*(10**(22))
    #number of time steps = number of output data length
    
    t_start = t_gps - seglen/2.
    t_end   = t_gps + time_shift - 1./srate
    adseglen = t_end - t_start
    t1 = np.linspace(t_start, t_end, int(adseglen*srate))
    #t_gps + time_shift = t_coalescence
    t = t_gps + time_shift - t1
    #shorten t array
    
    f = np.zeros(shape=srate*seglen)
    
    f[0:int(adseglen*srate)] = (256./5.)**(-3./8.)*1./(np.pi)*(G*chirpmass/(c**3.))**(-5./8.)*t**(-3./8.)
    #f = 20HZ at t1 = t_gps
    
    hp = np.zeros(shape=np.size(f))
    hc = np.zeros(shape=np.size(f))
    
    hp[0:int(adseglen*srate)] = 4/d*(G*chirpmass/(c**2))**(5/3)*(np.pi*f[0:int(adseglen*srate)]/c)**(2/3)*(np.cos(2*np.pi*f[0:int(adseglen*srate)]*t[:] + phiRef))*(np.cos(jota)**2+1.)/2
    hc[0:int(adseglen*srate)] = 4/d*(G*chirpmass/(c**2))**(5/3)*(np.pi*f[0:int(adseglen*srate)]/c)**(2/3)*(np.sin(2*np.pi*f[0:int(adseglen*srate)]*t[:] + phiRef))*np.cos(jota)
    
    
    #strain for both TT-Polarizations 
    #h computed in SI-Units
    #I lost a factor ~3 somehow...
    #f(t) have to be corrected!!! (will do on 30.05)
    return hp, hc