import numpy as np
from photon_journey.photon import Photon
def launch_a_photon(beam_radius,angle,z):
    #default to luanching photons into x-y plane at position z with an angle
    #stokes is defined as 1,1,0,0
    r=beam_radius*np.sqrt(-np.log(np.random.rand()))
    beta=np.pi*2*np.random.rand()
    x=r*np.cos(beta)
    y=r*np.sin(beta)

    dx= np.cos(angle)
    dy= 0
    dz= np.cos(angle)
    p= Photon(position=np.array([x,y,z]),direction=np.array([dx,dy,dz]),energy=1,stokes=np.array([1.0,1.0,0.0,0.0]))
    return p


