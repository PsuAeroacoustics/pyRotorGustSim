#!/usr/bin/env python3

from geometry import *
from bemt import *
import numpy as np
import aerosandbox as asb
import scipy.optimize as opt
import matplotlib.pyplot as plt
# import neuralfoil as nf
# from xfoil import XFoil
# from xfoil import model

#%%

def compute_aero(geom_params,input_params,res_param,observer_params,acs_params,saved_params,filt):
    
    Nb = geom_params['number_of_blades']
    R = geom_params['radius']
    e = geom_params['r_c']*R
    AR = geom_params['AR']
    TR = geom_params['TR']
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

    af = asb.Airfoil(geom_params['airfoil'])
    # af.coordinates = af.repanel(n_points_per_side = int(saved_params['airfoil_points']/2)).coordinates*saved_params['c']

    # aoa_polar = np.arange((20-(-5))/.5+1)*.5-5
    # Re = 0.75*R*omega*c/nu
    # M = 0.75*omega*R/sos
    # polar = af.get_aero_from_neuralfoil(alpha=aoa_polar, Re=Re,mach =M,model_size='xlarge')
    # aoa_cp_dist = 0.5*(aoa_polar[polar['CL'].argmax()]+aoa_polar[np.abs(polar['CL']).argmin()])

    # xf = XFoil()
    # xf.airfoil = model.Airfoil(x = af.coordinates[:,0],y = af.coordinates[:,1])    
    # xf.max_iter = 100
    # xf.Re = Re
    # xf.M = M
    # cl,cd = xf.a(aoa_cp_dist, as_dict=False)[:2]

    # cp = xf.get_cp_distribution()[-1]
    # cl = np.trapz(cp,x = af.coordinates[:,0],axis = -1)*np.cos(aoa_cp_dist*np.pi/180)+np.trapz(cp,x = af.coordinates[:,-1],axis = -1)*np.sin(aoa_cp_dist*np.pi/180)
    # cp_dist = cp/cl

    # out = xf.a(5, as_dict=True)
    # x,y,cp = xf.get_cp_distribution()
    # np.trapz(cp)
    # th_u = (np.arctan2(np.gradient(y),-np.gradient(x)))[:int(len(cp)/2)]
    # th_l = (np.arctan2(-np.gradient(y),np.gradient(x)))[int(len(cp)/2):]
    # # np.trapz(cp,x = x)*np.cos(5*np.pi/180)+np.trapz(cp,x = y)*np.sin(5*np.pi/180)
    # cn = np.trapz((cp[int(len(cp)/2)+1:]-cp[:int(len(cp)/2)][::-1]),x = x[:int(len(cp)/2)][::-1])
    # ca = np.trapz((cp[:int(len(cp)/2)][::-1]-cp[int(len(cp)/2)+1:]),x = y[int(len(cp)/2)+1:])

    # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
    # ax.plot(x,y)
    # ax.plot(x[:int(len(cp)/2)],cp[:int(len(cp)/2)])
    # ax.plot(x[int(len(cp)/2):],cp[int(len(cp)/2):])


    dpsi = input_params['computational_params']['d_psi']*np.pi/180
    iterations = int((input_params['computational_params']['number_of_revs']+1)*(2*np.pi)/dpsi)

    # n_revs = np.ceil(omega/(2*np.pi*res_param['df']))
    # iterations = int(n_revs*2*np.pi/dpsi*2)
    psi = np.arange(iterations)*dpsi
    
    dt = dpsi/omega
    t = np.arange(iterations)*dt
    

    # initialize aircraft and rotor object
    atmos = Atmosphere()
    ac = aircraft(1)
    ac.rotors = [rotor(Nb = Nb,R = R,e = e,AR = AR,TR = TR,th0 = th0,th_tw = th_tw,N_elements = N_elements,af = af,Cl_a=Cl_a,origin = origin,omega = omega,V_c = V_c,C_T_target = C_T_target,atmos=atmos) for r_iter in range(ac.N_rotor)]

    # trims rotor in regular hovering flight
    # if filt:
    #     xtr_upper = res_param['c_extents'][0]
    # else:
    #     xtr_upper = 1
    xtr_upper = res_param['c_extents'][0]

    th0 = opt.newton(trim,x0 = ac.rotors[0].th0,args=(ac.rotors[0],ac.rotors[0].blades[0],xtr_upper),tol=5e-6,full_output=False)
    ac.rotors[0].th0 = th0
    ac.rotors[0].blades[0].set_twist(ac.rotors[0])
    ac.rotors[0].blades[0].set_loads()
    lam_bemt = ac.rotors[0].blades[0].lam
    # print(ac.rotors[0].blades[0].c)
    
    #%%

    if input_params['computational_params']['unsteady_loading']:


        assert ( 'azimuthal_location' in input_params['gust_params'] or 'gust_end_pnts' in input_params['gust_params'] or 'r_trace' in input_params['gust_params']), "Please specify one of the following in the 'gust_params' entry of the parameter file: 'azimuthal_location', 'gust_end_pnts', 'r_trace'" 

        if 'azimuthal_location' in input_params['gust_params']:
            psi_g = input_params['gust_params']['azimuthal_location']*np.pi/180
            h = ((np.expand_dims(psi,axis = -1)%(2*np.pi)-psi_g)%(2*np.pi))*ac.rotors[0].blades[0].r

        elif 'gust_end_pnts' in input_params['gust_params']:
            gust_end_pnts = np.array(input_params['gust_params']['gust_end_pnts'])/R
            gust_psi_lim = np.round(np.sort([np.arctan2(gust_end_pnts[1],gust_end_pnts[0]),np.arctan2(gust_end_pnts[-1],gust_end_pnts[-2])])/dpsi).astype(int)
            gust_psi_lim_ind = slice(gust_psi_lim[0],gust_psi_lim[-1])

            gust_r_lim = np.sort([np.linalg.norm(gust_end_pnts[:2]),np.linalg.norm(gust_end_pnts[-2:])])
            gust_r_lim_ind = np.abs(np.expand_dims(gust_r_lim,axis = -1)-ac.rotors[0].blades[0].r).argmin(axis = -1)
            gust_r_lim_ind = slice(gust_r_lim_ind[0],gust_r_lim_ind[-1])

            x_b = np.expand_dims(ac.rotors[0].blades[0].r[gust_r_lim_ind],axis = -1)*np.cos(psi)[gust_psi_lim_ind]
            y_b = np.expand_dims(ac.rotors[0].blades[0].r[gust_r_lim_ind],axis = -1)*np.sin(psi)[gust_psi_lim_ind]

            m = np.diff(gust_end_pnts[1::2])/np.diff(gust_end_pnts[::2])
            y_g = m*(x_b-gust_end_pnts[0])+gust_end_pnts[1]
            gust_ind = np.abs(y_b-y_g).argmin(axis = 1)

            x_g = ac.rotors[0].blades[0].r[gust_r_lim_ind]*np.cos(psi)[gust_psi_lim_ind][gust_ind]
            y_g = ac.rotors[0].blades[0].r[gust_r_lim_ind]*np.sin(psi)[gust_psi_lim_ind][gust_ind]

            psi_g = np.arctan2(y_g,x_g)

            h = np.zeros((iterations,N_elements))
            h[:,gust_r_lim_ind] = ((np.expand_dims(psi,axis = -1)%(2*np.pi)-psi_g)%(2*np.pi))*ac.rotors[0].blades[0].r[gust_r_lim_ind]

        elif 'r_trace' in input_params['gust_params']:
            
            psi_g = np.arctan2(1,input_params['gust_params']['r_trace'])
            x0 = input_params['gust_params']['r_trace']/np.cos(psi_g)
            psi_lim = slice(int(np.arctan2(np.sqrt(1-(x0/2-x0)**2),-x0/2)/dpsi),int(np.pi/dpsi))
            x_g = np.cos(psi[psi_lim][::-1])+x0
            y_g = np.sin(psi[psi_lim][::-1])

            r_g = np.linalg.norm((x_g,y_g),axis = 0)
            psi_g = np.arctan2(y_g,x_g)

            r_ind = np.abs(np.expand_dims(ac.rotors[0].blades[0].r,axis = -1)-r_g).argmin(axis = -1)
            psi_g =psi[np.abs(np.expand_dims(psi,axis = -1) -psi_g[r_ind]).argmin(axis = 0)]
            x_g = ac.rotors[0].blades[0].r*np.cos(psi_g)
            y_g = ac.rotors[0].blades[0].r*np.sin(psi_g)

            h = ((np.expand_dims(psi,axis = -1)%(2*np.pi)-psi_g)%(2*np.pi))*ac.rotors[0].blades[0].r

        #$$

        gamma_v = input_params['gust_params']['strength']
        y = -input_params['gust_params']['peak_location']/AR
        n = 2
        r_c = 0.05/AR

        # Nb = 4
        # sigma = Nb*c/(np.pi*R) 
        # CT = 0.012
        # gamma_v = 2*CT/sigma
        # y = -0.25*c/R
        # h = np.abs(y)
        # h = (np.arange(200)*(2/AR)/(200-1)-(2/AR)/2)
        # x = h
        r = np.sqrt(h**2+y**2)
        lam_gust= ((R/AR)/R*gamma_v/(2*np.pi)*(r/(r_c**(2*n)+r**(2*n))**(1/n)))*h/r

        # lam_gust = np.zeros(lam_gust_1.shape)
        # ind = np.abs(0.75-ac.rotors[0].blades[0].r).argmin()
        # lam_gust[:,ind] = 15*lam_gust_1[:,ind]

        v_gust = lam_gust*omega*R
        ind = np.where(v_gust==0)
        v_gust[ind] = v_gust[tuple((ind[0]-1,ind[1]))]
        # v_gust_2 = 2*omega*R*c*CT/sigma*1/(2*np.pi)*(1/((r_c*R)**(2*n)+(r*R)**(2*n)**(1/n))*h*R)
        # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
        # ax.plot(h[:int(2*np.pi/dpsi),-1]*AR,v_gust[:int(2*np.pi/dpsi),-1],label = 'gust velocity')
        # # ax.set_xlim([-0.2,1.4])
        # # ax.set_xticks(np.arange(9)*0.2-0.2)
        # # ax.set_ylim([0,160])
        # ax.grid()


        # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
        # ax.plot(psi[:int(2*np.pi/dpsi)],v_gust[:int(2*np.pi/dpsi)]*3.28084,label = 'gust velocity')
        # ax.legend(np.round(ac.rotors[0].blades[0].r[::4],2))
        # # ax.set_xlim([-0.2,1.4])
        # # ax.set_ylim([0,160])
        # ax.grid()

        # fig,ax = plt.subplots(subplot_kw=dict(projection = 'polar'))
        # levels = np.linspace(0,35,21)
        # dist = ax.contourf(psi[:int(2*np.pi/dpsi)],ac.rotors[0].blades[0].r,v_gust[:int(2*np.pi/dpsi)].T,levels = levels)
        # cbar = fig.colorbar(dist,pad = .1)
        # cbar.ax.set_ylabel(r'$\alpha \ [deg]$')


         # Indicial response function coefficients and exponents (these are derived from CFD data and given by Leishman)
        A1 = 0.67
        b1 = .1753
        A2 = 0.33
        b2 = 1.637

        # total inflow ratio accounting for the gust contributions
        lam = lam_gust+lam_bemt
        U = np.sqrt(lam**2+ac.rotors[0].blades[0].r**2)
        beta = np.sqrt(1-(U*omega*R/sos)**2)
        
        # non-dimensionalized distance in terms of half chords
        s = omega*R*ac.rotors[0].blades[0].r*np.expand_dims(t,axis = -1)/(ac.rotors[0].blades[0].c/2)
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
        # aoa_eff = aoa_eff+np.expand_dims(ac.rotors[0].blades[0].aoa,axis = 0)
        # aoa_eff = (ac.rotors[0].blades[0].th - lam_gust/ac.rotors[0].blades[0].r)
        # new effective inflow angle
        phi_eff = ac.rotors[0].blades[0].th - aoa_eff
        # spanwise Reynold's and Mach numbers
        Re = np.ones(aoa_eff.shape)*ac.rotors[0].blades[0].Re
        M = np.ones(aoa_eff.shape)*ac.rotors[0].blades[0].M
        CL,CD = get_af_coeffs(ac.rotors[0].blades[0].af,aoa_eff*180/np.pi,Re,M)

        # sectional axial and normal force coefficients
        dCz = CL*np.cos(phi_eff)-CD*np.sin(phi_eff)
        dCx = CL*np.sin(phi_eff)+CD*np.cos(phi_eff)

        
        # cl = np.zeros(360)
        # cd = np.zeros(360)
        # cp = np.zeros((360,len(af.coordinates)))

        # for i in range(360):
        #     xf.Re = Re[i,int(0.75*N_elements)] 
        #     xf.M = M[i,int(0.75*N_elements)]
        #     cl[i],cd[i] = xf.a(aoa_eff[i,int(0.75*N_elements)]*180/np.pi, as_dict=False)[:2]
        #     cp[i] = xf.get_cp_distribution()[-1]

        # cl_cp = np.trapz(cp,x = af.coordinates[:,0],axis = -1)*np.cos(aoa_eff[:360,int(0.75*N_elements)])+np.trapz(cp,x = af.coordinates[:,-1],axis = -1)*np.sin(aoa_eff[:360,int(0.75*N_elements)])

        # cp = cp_dist*np.expand_dims(CL,axis = -1)
        # cl_cp = np.trapz(cp,x = af.coordinates[:,0],axis = -1)*np.cos(aoa_cp_dist*np.pi/180)+np.trapz(cp,x = af.coordinates[:,-1],axis = -1)*np.sin(aoa_cp_dist*np.pi/180)
        # dCz = cl_cp*np.cos(phi_eff)-CD*np.sin(phi_eff)


        # dCT_cp = np.expand_dims(0.5*ac.rotors[0].c/(np.pi*ac.rotors[0].R)*U**2,axis = -1)*cp
        # dFz_cp = dCT_cp*atmos.rho*np.pi*ac.rotors[0].R*(ac.rotors[0].omega*ac.rotors[0].R)**2
        # dFz_cp = np.trapz(dFz_cp,x = af.coordinates[:,0],axis = -1)*np.cos(aoa_cp_dist*np.pi/180)+np.trapz(dFz_cp,x = af.coordinates[:,-1],axis = -1)*np.sin(aoa_cp_dist*np.pi/180)


        # ax.plot(cl_cp[:360,int(0.75*N_elements)])
        # ax.legend(['neuralfoil','xfoil','from cp'])

        # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
        # ax.plot(dFz[:360,int(0.75*N_elements)])
        # ax.plot(dFz_cp[:360,int(0.75*N_elements)])
        # ax.legend(['neuralfoil','xfoil','from cp'])


    else:
        U = ac.rotors[0].blades[0].U*np.ones((iterations,N_elements))
        dCz =(ac.rotors[0].blades[0].CL*np.cos(ac.rotors[0].blades[0].phi)-ac.rotors[0].blades[0].CD*np.sin(ac.rotors[0].blades[0].phi))*np.ones((iterations,N_elements))
        dCx = (ac.rotors[0].blades[0].CL*np.sin(ac.rotors[0].blades[0].phi)+ac.rotors[0].blades[0].CD*np.cos(ac.rotors[0].blades[0].phi))*np.ones((iterations,N_elements))

    # sectional thrust and power coefficients
    dCT = 0.5*ac.rotors[0].c/(np.pi*ac.rotors[0].R)*U**2*dCz
    dCP = 0.5*ac.rotors[0].c/(np.pi*ac.rotors[0].R)*ac.rotors[0].blades[0].r*U**2*dCx

    # sectional dimensionalized blade loads
    dFx= dCP*rho*np.pi*ac.rotors[0].R**2*(ac.rotors[0].omega*ac.rotors[0].R)**2
    dFz = dCT*rho*np.pi*ac.rotors[0].R*(ac.rotors[0].omega*ac.rotors[0].R)**2
    dFy = np.zeros(dFz.shape)
    loads = np.array([dFy,-dFx,dFz]).transpose(1,2,0)

    # lifting line nodes and normals
    lifting_line_nodes = np.expand_dims(np.array([ac.rotors[0].R*ac.rotors[0].blades[0].r,np.zeros(N_elements),np.zeros(N_elements)]).T,axis = 0)
    lifting_line_norms = np.expand_dims(np.array([np.zeros(N_elements),np.zeros(N_elements),np.ones(N_elements)]).T,axis = 0)

    if input_params['computational_params']['unsteady_loading']:
        saved_params.update({'t':t,'iterations':iterations,'r_elem':ac.rotors[0].blades[0].elems,'r':ac.rotors[0].blades[0].r,"th":ac.rotors[0].blades[0].th,'airfoil':geom_params['airfoil'],'airfoil_points':input_params['computational_params']['airfoil_elements'],'th_tw':th_tw,'TR':TR,'AR':AR,'R':R,'e':e,'c':ac.rotors[0].blades[0].c,'dpsi':dpsi,'dt':dt,'sos':sos,'N_elements':N_elements,'omega':omega,'psi':psi,'h_gust':h,'v_gust':v_gust,'U':U,'s':s,'th0':th0,'CL':CL,'CD':CD,'aoa':aoa_eff,'phi_eff':phi_eff,'dCT':dCT,'dCP':dCP,'loads':loads,'lifting_line_nodes':lifting_line_nodes,'lifting_line_norms':lifting_line_norms})
    else:
        saved_params.update({'t':t,'iterations':iterations,'r_elem':ac.rotors[0].blades[0].elems,'r':ac.rotors[0].blades[0].r,"th":ac.rotors[0].blades[0].th,'airfoil':geom_params['airfoil'],'airfoil_points':input_params['computational_params']['airfoil_elements'],'th_tw':th_tw,'TR':TR,'AR':AR,'R':R,'e':e,'c':ac.rotors[0].blades[0].c,'dpsi':dpsi,'dt':dt,'sos':sos,'N_elements':N_elements,'omega':omega,'psi':psi,'th0':th0,'CL':ac.rotors[0].blades[0].CL,'CD':ac.rotors[0].blades[0].CD,'aoa':ac.rotors[0].blades[0].aoa,'dCT':dCT,'dCP':dCP,'loads':loads,'lifting_line_nodes':lifting_line_nodes,'lifting_line_norms':lifting_line_norms})


