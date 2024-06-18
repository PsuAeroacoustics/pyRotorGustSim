#!/usr/bin/env python3

from geometry import *
from bemt import *
import numpy as np
import aerosandbox as asb
import scipy.optimize as opt


#%%

def compute_aero(geom_params,input_params,observer_params,acs_params,saved_params):
    
    Nb = geom_params['number_of_blades']
    R = geom_params['radius']
    e = geom_params['r_c']*R
    c = R/geom_params['AR']
    th_tw = geom_params['theta_tw']*np.pi/180
    th0 = geom_params['theta_initial']*np.pi/180
    Cl_a = 2*np.pi
    origin = [0,0,0]

    rho = input_params['flight_params']['density']
    sos = input_params['flight_params']['sos']
    nu = input_params['flight_params']['kinematic_viscosity']

    N_elements = input_params['computational_params']['spanwise_elements']
    C_T_target = input_params['flight_params']['C_T_target']
    omega = input_params['flight_params']['omega']
    V_c = 0
    af = asb.Airfoil("naca0009")
    dpsi = input_params['computational_params']['d_psi']*np.pi/180
    iterations = int(input_params['computational_params']['number_of_revs']*(2*np.pi)/dpsi)
    psi = np.arange(iterations)*dpsi
    dt = dpsi/omega
    t = np.arange(iterations)*dt
    
    gamma_vf = input_params['gust_params']['strength']/R
    r_vf = input_params['gust_params']['core_size']*c/R
    n = 2
    psi_gust = input_params['gust_params']['azimuthal_location']*np.pi/180

    # Indicial response function coefficients and exponents (these are derived from CFD data and given by Leishman)
    A1 = 0.67
    b1 = .1753
    A2 = 0.33
    b2 = 1.637

    # initialize aircraft and rotor object
    atmos = Atmosphere()
    ac = aircraft(1)
    ac.rotors = [rotor(Nb = Nb,R = R,e = e,c = c,th0 = th0,th_tw = th_tw,N_elements = N_elements,af = af,Cl_a=Cl_a,origin = origin,omega = omega,V_c = V_c,C_T_target = C_T_target,atmos=atmos) for r_iter in range(ac.N_rotor)]

    # trims rotor in regular hovering flight
    th0 = opt.newton(trim,x0 = ac.rotors[0].th0,args=(ac.rotors[0],ac.rotors[0].blades[0]),tol=5e-6,full_output=False)
    ac.rotors[0].th0 = th0
    ac.rotors[0].blades[0].set_twist(ac.rotors[0])
    ac.rotors[0].blades[0].set_loads()
    lam_bemt = ac.rotors[0].blades[0].lam

    # gust profile and induced velocity
    h = np.expand_dims(((psi%(2*np.pi)-psi_gust)%(2*np.pi)),axis = -1)*ac.rotors[0].blades[0].r
    v_gust = gamma_vf/(2*np.pi)*(h/(r_vf**(2*n)+(h)**(2*n))**(1/n))
    lam_gust = v_gust/(omega*R)
    
    # total inflow ratio accounting for the gust contributions
    lam = lam_bemt+lam_gust
    U = np.sqrt(lam**2+ac.rotors[0].blades[0].r**2)
    beta = np.sqrt(1-(U*omega*R/sos)**2)
    
    # non-dimensionalized distance in terms of half chords
    s = omega*R*ac.rotors[0].blades[0].r*np.expand_dims(t,axis = -1)/(c/2)
    ds = np.diff(s,axis = 0)[0]

    aoa_eff = np.zeros((iterations,N_elements))
    X_temp = np.zeros((N_elements))
    Y_temp = np.zeros((N_elements))

    for i in range(iterations):

        X = X_temp*np.exp(-b1*beta[i]**2*ds)+A1*omega*R*(lam_gust[i]-lam_gust[i-1])*np.exp(-b1*beta[i]**2*ds)**(1/2)
        Y = Y_temp*np.exp(-b2*beta[i]**2*ds)+A2*omega*R*(lam_gust[i]-lam_gust[i-1])*np.exp(-b2*beta[i]**2*ds)**(1/2)
        dCL = 2*np.pi/(beta[i]*U[i]*omega*R)*(lam_gust[i]*omega*R-X-Y)
        aoa_eff[i] = dCL/Cl_a
        X_temp = X
        Y_temp = Y

    # new effective inflow angle
    phi_eff = ac.rotors[0].blades[0].th - aoa_eff

    # spanwise Reynold's and Mach numbers
    Re = np.ones(aoa_eff.shape)*ac.rotors[0].blades[0].Re
    M = np.ones(aoa_eff.shape)*ac.rotors[0].blades[0].M

    # computes sectional airfoil coefficients
    CL,CD = get_af_coeffs(ac.rotors[0].blades[0].af,aoa_eff*180/np.pi,Re,M)

    # sectional axial and normal force coefficients
    dCz = CL*np.cos(phi_eff)-CD*np.sin(phi_eff)
    dCx = CL*np.sin(phi_eff)+CD*np.cos(phi_eff)

    # sectional thrust and power coefficients
    dCT = 0.5*ac.rotors[0].c/(np.pi*ac.rotors[0].R)*U**2*dCz
    dCP = 0.5*ac.rotors[0].c/(np.pi*ac.rotors[0].R)*ac.rotors[0].blades[0].r*U**2*dCx

    # sectional dimensionalized blade loads
    dFx= dCP*atmos.rho*np.pi*ac.rotors[0].R**2*(ac.rotors[0].omega*ac.rotors[0].R)**2
    dFz = dCT*atmos.rho*np.pi*ac.rotors[0].R*(ac.rotors[0].omega*ac.rotors[0].R)**2
    dFy = np.zeros(dFz.shape)
    loads = np.array([dFy,-dFx,dFz]).transpose(1,2,0)

    # lifting line nodes and normals
    lifting_line_nodes = np.expand_dims(np.array([ac.rotors[0].R*ac.rotors[0].blades[0].r,np.zeros(N_elements),np.zeros(N_elements)]).T,axis = 0)
    lifting_line_norms = np.expand_dims(np.array([np.zeros(N_elements),np.zeros(N_elements),np.ones(N_elements)]).T,axis = 0)

    saved_params.update({'t':t,'iterations':iterations,'omega':omega,'psi':psi,'h_gust':h,'v_gust':v_gust,'s':s,'th0':th0,'CL':CL,'CD':CD,'aoa_eff':aoa_eff,'dCT':dCT,'dCP':dCP,'loads':loads,'lifting_line_nodes':lifting_line_nodes,'lifting_line_norms':lifting_line_norms})


