#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
# import neuralfoil as nf
import aerosandbox as asb
import scipy.optimize as opt

#%%

rho = 1.125
sos = 340
nu = 14.88e-6

#%%

Nb = 4
R = 2
e = .22*R
c = R/12.5
th_tw = -0*np.pi/180
th0 = 10*np.pi/180
N_elements = 48
Cl_a = 2*np.pi

# parameterizes airfoil for neuralfoil
airfoil = asb.Airfoil("naca0012")

#%%

omega = 109.12
V_c = 0
C_T_target = 0.008

#%%

r_elem = (np.arange(N_elements+1)*(R-e)/(N_elements)+e)/R
r = 0.5*(r_elem[1:]+r_elem[:-1])

th = th0+r*th_tw

sigma = Nb*c/(np.pi*R)

lam_c = V_c/(omega*R)

#%%

Re = omega*R*r*c/nu
M = omega*R*r/sos

#%%

def get_loads(lam,th):

    U = np.sqrt(lam**2+r**2)
    phi = np.arctan2(lam,r)
    aoa = th-phi

    aero = airfoil.get_aero_from_neuralfoil(alpha=aoa*180/np.pi, Re=Re,mach = M,model_size='large')

    dCz = aero['CL']*np.cos(phi)-aero['CD']*np.sin(phi)
    dCx = aero['CL']*np.sin(phi)+aero['CD']*np.cos(phi)

    dCT = sigma/2*U**2*dCz
    dCP = sigma/2*r*U**2*dCx

    CT = np.trapz(dCT,x = r)    
    CP = np.trapz(dCP,x = r)    

    return aoa, dCT, dCP, CT, CP


def solve_lam(lam_init,th):

    phi = np.arctan2(lam_init,r)
    aoa = th-phi

    f = Nb/2*(1-r)/(r*phi)
    F = (2/np.pi)*np.arccos(np.exp(-f))

    aero = airfoil.get_aero_from_neuralfoil(alpha=aoa*180/np.pi, Re=Re,mach = M,model_size='large')
    dCz = aero['CL']*np.cos(phi)-aero['CD']*np.sin(phi)

    lam = (4*lam_c*r*F+np.sqrt((4*lam_c*r*F)**2-4*(4*r*F-sigma/2*dCz)*(-sigma/2*dCz*r**2)))/(2*(4*r*F-sigma/2*dCz))

    err = (lam-lam_init)/lam
    err[np.isnan(err)] = 1

    return err

def solve_lam2(lam_init,th):

    phi = np.arctan2(lam_init,r)
    f = Nb/2*(1-r)/(r*phi)
    F = (2/np.pi)*np.arccos(np.exp(-f))
    lam = sigma*Cl_a/(16*F)*(np.sqrt(1+32*F/(sigma*Cl_a)*th*r)-1)

    return lam



def trim(th0):

    th = th0+r*th_tw

    lam_init = np.sqrt((sigma*Cl_a/16-lam_c/2)**2+sigma*Cl_a/8*th*r)-(sigma*Cl_a/16-lam_c/2)
    lam = opt.newton(solve_lam,x0 = lam_init,args=(th,))

    aoa, dCT, dCP, CT, CP  = get_loads(lam,th)

    err = abs(C_T_target-CT)/C_T_target

    if np.any(np.isnan(lam)):
        err = 1

    print(err)

    return err


def trim2(th0):
    th = th0+r*th_tw

    lam_init = np.sqrt((sigma*Cl_a/16-lam_c/2)**2+sigma*Cl_a/8*th*r)-(sigma*Cl_a/16-lam_c/2)
    lam = opt.fixed_point(solve_lam2,x0 = lam_init,args=(th,))
    # lam = sigma*Cl_a/16*(np.sqrt(1+32/(sigma*Cl_a)*th*r)-1)
    aoa, dCT, dCP, CT, CP  = get_loads(lam,th)
    err = abs(C_T_target-CT)/C_T_target
    return(err)

sol = opt.newton(trim,x0 = th0,tol=5e-5,full_output=True)

sol2 = opt.newton(trim2,x0 = th0,tol=5e-5,full_output=True)

th = sol[0]+r*th_tw
lam_init = np.sqrt((sigma*Cl_a/16-lam_c/2)**2+sigma*Cl_a/8*th*r)-(sigma*Cl_a/16-lam_c/2)
lam = opt.newton(solve_lam,x0 = lam_init,args=(th,),full_output=False)
aoa, dCT, dCP, CT, CP = get_loads(lam,th)

th2 = sol2[0]+r*th_tw
lam_init = sigma*Cl_a/16*(np.sqrt(1+32/(sigma*Cl_a)*th2*r)-1)
lam2 = opt.fixed_point(solve_lam2,x0 = lam_init,args=(th2,))
aoa, dCT2, dCP2, CT2, CP2 = get_loads(lam2,th2)


#%%

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
ax.plot(r,lam)
ax.plot(r,lam2)

ax.set_xlabel('r/R')
ax.grid()
ax.set_ylabel('$\lambda$')
ax.set_xlim([0,1])
ax.set_ylim([0,.1])


fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
ax.plot(r,dCT)
ax.plot(r,dCT2)

ax.set_xlabel('r/R')
ax.grid()
# ax.set_ylabel('$\dC_T$')
ax.set_xlim([0,1])
# ax.set_ylim([0,.03])
