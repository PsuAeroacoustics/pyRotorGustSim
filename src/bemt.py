import numpy as np
import scipy.optimize as opt

#%%

def solve_lam(lam_init,rotor, blade):
    blade.phi = np.arctan2(lam_init,blade.r)
    blade.aoa = blade.th-blade.phi
    f =rotor.Nb/2*(1-blade.r)/(blade.r*blade.phi)
    F = (2/np.pi)*np.arccos(np.exp(-f))
    lam = rotor.sigma*blade.Cl_a/(16*F)*(np.sqrt(1+32*F/(rotor.sigma*blade.Cl_a)*blade.th*blade.r)-1)
    return lam

    
def trim(th0,rotor,blade):
    
    rotor.th0 = th0
    blade.th = rotor.th0+blade.r*rotor.th_tw

    lam_init = np.sqrt((rotor.sigma*blade.Cl_a/16-rotor.lam_c/2)**2+rotor.sigma*blade.Cl_a/8*blade.th*blade.r)-(rotor.sigma*blade.Cl_a/16-rotor.lam_c/2)
    blade.lam = opt.fixed_point(solve_lam,x0 = lam_init,args=(rotor,blade))

    # blade.lam = opt.newton(solve_lam,x0 = blade.lam,args=(rotor,blade))

    blade.set_loads()
    rotor.CT = rotor.Nb*np.trapz(blade.dCT,x = blade.r)

    err = abs(rotor.C_T_target-rotor.CT)/rotor.C_T_target

    if np.any(np.isnan(blade.lam)):
        err = 1

    print(err)
    return err