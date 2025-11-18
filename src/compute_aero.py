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

    dpsi = input_params['computational_params']['d_psi']*np.pi/180
    iterations = int((input_params['computational_params']['number_of_revs']+1)*(2*np.pi)/dpsi)

    psi = np.arange(iterations)*dpsi
    
    dt = dpsi/omega
    t = np.arange(iterations)*dt
    

    # initialize aircraft and rotor object
    atmos = Atmosphere()
    ac = aircraft(1)
    ac.rotors = [rotor(Nb = Nb,R = R,e = e,AR = AR,TR = TR,th0 = th0,th_tw = th_tw,N_elements = N_elements,af = af,Cl_a=Cl_a,origin = origin,omega = omega,V_c = V_c,C_T_target = C_T_target,atmos=atmos) for r_iter in range(ac.N_rotor)]

    xtr_upper = res_param['c_extents'][0]

    th0 = opt.newton(trim,x0 = ac.rotors[0].th0,args=(ac.rotors[0],ac.rotors[0].blades[0],xtr_upper),tol=5e-6,full_output=False)
    ac.rotors[0].th0 = th0
    ac.rotors[0].blades[0].set_twist(ac.rotors[0])
    ac.rotors[0].blades[0].set_loads()
    lam_bemt = ac.rotors[0].blades[0].lam
    # print(ac.rotors[0].blades[0].c)
    
    #%%

    if input_params['computational_params']['unsteady_loading']:


        #$$
        extents_ind = np.abs(ac.rotors[0].blades[0].r[:,None]-input_params['af_params']['extents']).argmin(0)
        D = input_params['af_params']['radius']/R
        y = (input_params['af_params']['seperation_distance']/R+np.sign(input_params['af_params']['seperation_distance'])*D)
        x = np.roll((psi[:,None]%(np.pi*2)-np.pi),-input_params['af_params']['azimuthal_location']*np.pi/180/dpsi,axis = 0)*ac.rotors[0].blades[0].r[extents_ind[0]:extents_ind[1]]
        # r = np.sqrt((x**2+y**2))
        # th = np.arctan2(-y,-x)+np.pi/2

        # lam_r = lam_bemt*np.cos(th)*(1-(D/r)**2)
        # lam_th = -lam_bemt*np.sin(th)*(1+(D/r)**2)
        # lam_x = lam_r*np.sin(th)+lam_th*np.cos(th)
        # lam_y = lam_r*np.cos(th)-lam_th*np.sin(th)

        z_v =x+1j*y
        w = -1j*lam_bemt[extents_ind[0]:extents_ind[1]]*(1+(D/z_v)**2)
        lam_x,lam_y = np.real(w), -np.imag(w)

        # lam_af = lam_y-lam_bemt
        # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
        # plt.plot(psi, lam_y[:,int(N_elements*0.7)])
        # # plt.plot(psi, np.sqrt(lam_y**2+lam_x**2)[:,int(N_elements*0.7)],linestyle = '-.')
        # # plt.plot(psi, np.sqrt(V_x**2+V_y**2)[:,int(N_elements*0.7)],linestyle = '-.')
        # plt.plot(psi, V_y[:,int(N_elements*0.7)],linestyle = '-.')

        # V = 1
        # D = 2
        # x = np.arange(100)/99*12-6
        # y = np.arange(120)/119*14-7
        # x_mesh,y_mesh = np.meshgrid(x,y)
        # r = np.sqrt(x_mesh**2+y_mesh**2)
        # th = np.arctan(y_mesh/x_mesh)
     
        # V_r = V*np.cos(th)*(1-(D/r)**2)
        # V_th = -V*np.sin(th)*(1+(D/r)**2)
        # V_y = V_r*np.sin(th)+V_th*np.cos(th)
        # V_x = V_r*np.cos(th)-V_th*np.sin(th)
    
        # levels = np.linspace(-5,5,100)
        # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
        # # plt.contourf(x_mesh, y_mesh, np.sqrt(V_y**2+V_x**2),cmap = 'inferno',levels = levels)
        # plt.contourf(x_mesh, y_mesh, V_x,cmap = 'inferno',levels = levels)

        # z = -y[:,None]+1j*x
        # w = V*(1-(D/z)**2)
        # V_x,V_y = np.real(w), -np.imag(w).T

        # levels = np.linspace(-5,5,100)
        # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
        # # plt.contourf(x_mesh, y_mesh, np.sqrt(V_y**2+V_x**2),cmap = 'inferno',levels = levels)
        # # plt.contourf(x_mesh, y_mesh, np.sqrt(np.real(w)**2+np.imag(w)**2).T,cmap = 'inferno',levels = levels)
        # plt.contourf(x_mesh, y_mesh, np.real(w),cmap = 'inferno',levels = levels)
        # ax.set_ylim([-6,6])



         # Indicial response function coefficients and exponents (these are derived from CFD data and given by Leishman)
        A1 = 0.67
        b1 = .1753
        A2 = 0.33
        b2 = 1.637

        #%%
        # A1 = 0.5
        # b1 = .13
        # A2 = 0.5
        # b2 = 1.0

        
        # total inflow ratio accounting for the gust contributions
        aoa = (ac.rotors[0].blades[0].th[extents_ind[0]:extents_ind[1]]-np.arctan2(lam_y,ac.rotors[0].blades[0].r[extents_ind[0]:extents_ind[1]]+lam_x))
        U = np.sqrt(lam_y**2+(ac.rotors[0].blades[0].r[extents_ind[0]:extents_ind[1]]+lam_x)**2)
        beta = np.sqrt(1-(U*omega*R/sos)**2)
        
        # non-dimensionalized distance in terms of half chords
        s = omega*R*(ac.rotors[0].blades[0].r/(ac.rotors[0].blades[0].c/2))[extents_ind[0]:extents_ind[1]]*t[:,None]
        ds = np.diff(s[:2],axis = 0)[0]

        aoa_eff = (ac.rotors[0].blades[0].aoa[:,None]*np.ones(iterations)).T
        X_temp = np.zeros((extents_ind[1]-extents_ind[0]))
        Y_temp = np.zeros((extents_ind[1]-extents_ind[0]))

        for i in range(iterations):

            X = X_temp*np.exp(-b1*ds)+A1*(aoa[i]-aoa[i-1])*np.exp(-b1*2*ds/2)
            Y = Y_temp*np.exp(-b2*ds)+A2*(aoa[i]-aoa[i-1])*np.exp(-b2*2*ds/2)
            aoa_eff[i,extents_ind[0]:extents_ind[1]] = aoa[i]-X-Y
            X_temp = X
            Y_temp = Y

            # X = X_temp*np.exp(-b1*beta[i]**2*ds)+A1*omega*R*(lam_y[i,extents_ind[0]:extents_ind[1]]-lam_y[i-1,extents_ind[0]:extents_ind[1]])*np.exp(-b1*beta[i]**2*ds)**(1/2)
            # Y = Y_temp*np.exp(-b2*beta[i]**2*ds)+A2*omega*R*(lam_y[i,extents_ind[0]:extents_ind[1]]-lam_y[i-1,extents_ind[0]:extents_ind[1]])*np.exp(-b2*beta[i]**2*ds)**(1/2)
            # aoa_eff[i,extents_ind[0]:extents_ind[1]] = 1/(beta[i]*U[i]*omega*R)*(lam_y[i,extents_ind[0]:extents_ind[1]]*omega*R-X-Y)
            # X_temp = X
            # Y_temp = Y

        # new effective inflow angle
        phi_eff = ac.rotors[0].blades[0].th - aoa_eff
        lam_eff = (ac.rotors[0].blades[0].r*np.tan(phi_eff))
        U = np.sqrt(lam_eff**2+ac.rotors[0].blades[0].r**2)
        lam_eff = lam_eff[:,extents_ind[0]:extents_ind[1]]

        # spanwise Reynold's and Mach numbers
        Re = np.ones(aoa_eff.shape)*ac.rotors[0].blades[0].Re
        M = np.ones(aoa_eff.shape)*ac.rotors[0].blades[0].M
        CL,CD = get_af_coeffs(ac.rotors[0].blades[0].af,aoa_eff*180/np.pi,Re,M)
        # sectional axial and normal force coefficients
        dCz = CL*np.cos(phi_eff)-CD*np.sin(phi_eff)
        dCx = CL*np.sin(phi_eff)+CD*np.cos(phi_eff)

        #%%

        gamma = (0.5*ac.rotors[0].blades[0].c/ac.rotors[0].blades[0].R*U*CL)[:,extents_ind[0]:extents_ind[1]]
        th = np.arange(100)/100*2*np.pi
        dth = np.diff(th[:2])
        z_af = D*np.exp(1j*th)[:,None,None]


        w_af = -1j*lam_eff*(1+(D/z_af)**2)-1j*gamma/(2*np.pi)*(1/(z_af-z_v)-1/(z_af-D**2/np.conj(z_v))+1/z_af)
        # w_af = -1j*lam_eff*(1+(D/z_af)**2)-1j*gamma/(2*np.pi)*(1/(z_af-z_v)-1/(z_af-D**2/np.conj(z_v)))

        lam_x_af,lam_y_af = np.real(w_af), -np.imag(w_af)
        lam_r_af, lam_th_af =np.real(w_af.T*np.exp(1j*th)).T,-np.imag(w_af.T*np.exp(1j*th)).T

        dth_dt = np.real(-1j*np.gradient(lam_eff,edge_order=2,axis = 0)/dt*(z_af-D**2/z_af)
                         -(1j*gamma/(2*np.pi)*(-np.gradient(z_v,edge_order=2,axis = 0)/dt/(z_af-z_v)
                                              +np.gradient(D**2/np.conj(z_v),edge_order=2,axis = 0)/dt/(z_af-D**2/np.conj(z_v))
                                              -np.gradient(np.conj(z_v),edge_order=2,axis = 0)/dt/np.conj(z_v))
                                +1j*np.gradient(gamma,edge_order=2,axis = 0)/dt/(2*np.pi)*(np.log(z_af-z_v)-np.log(z_af-D**2/np.conj(z_v))-np.log(-np.conj(z_v))+np.log(z_af))))

        # dth_dt = np.real(-1j*np.gradient(lam_eff,edge_order=2,axis = 0)/dt*(z_af-D**2/z_af)
        #                  -(1j*gamma/(2*np.pi)*(-np.gradient(z_v,edge_order=2,axis = 0)/dt/(z_af-z_v)
        #                                       +np.gradient(D**2/np.conj(z_v),edge_order=2,axis = 0)/dt/(1-D**2/np.conj(z_v)))
        #                         +1j*np.gradient(gamma,edge_order=2,axis = 0)/dt/(2*np.pi)*(np.log(z_af-z_v)+np.log(z_af-D**2/np.conj(z_v)))))

        # p = 0.5*rho*(lam_y**2-abs(w_cyl)**2)-rho*np.real(dth_dt)
        # dth_dt = np.real(-1j*np.gradient(lam_eff,edge_order=2,axis = 0)/dt*(z_af-D**2/z_af)
        #             -1j*gamma/(2*np.pi)*(-np.gradient(z_v,edge_order=2,axis = 0)/dt/(z_af-z_v)
        #                                 +np.gradient(D**2/np.conj(z_v),edge_order=2,axis = 0)/dt/(1-D**2/np.conj(z_v))
        #                                 ))

        cp = (1-(abs(w_af)/lam_eff)**2)-2/lam_eff**2*dth_dt
        neg_cp_ind = np.where(cp<0)
        neg_cp_ind_mod = (neg_cp_ind[0],neg_cp_ind[1]-2,neg_cp_ind[2])
        cp[neg_cp_ind] = cp[neg_cp_ind_mod]

        dFz_af = 0.5*rho*(lam_eff*omega*R)**2*D*R**2*np.trapezoid(cp.T*np.sin(th),th,axis = -1).T*np.diff(ac.rotors[0].blades[0].r[:2])
        dFx_af = 0.5*rho*(lam_eff*omega*R)**2*D*R**2*np.trapezoid(cp.T*np.cos(th),th,axis = -1).T*np.diff(ac.rotors[0].blades[0].r[:2])

        # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))        
        # # ax.plot(psi,lam_eff[:,-2])
        # ax.plot(psi,cp[::25,:,-2].T)

        # #%%
        # CD_af = 1.2
        # dFz_af = 200*0.5*rho*(lam_eff*omega*R)**2*CD_af*2*input_params['af_params']['radius']
            # np.trapezoid(dFz_af,x = ac.rotors[0].blades[0].r*ac.rotors[0].blades[0].R,axis = -1).max()
        # lifting line nodes and normals
        af_nodes = np.expand_dims(np.array([ac.rotors[0].R*ac.rotors[0].blades[0].r[extents_ind[0]:extents_ind[1]]*np.cos(input_params['af_params']['azimuthal_location']*np.pi/180),ac.rotors[0].R*ac.rotors[0].blades[0].r[extents_ind[0]:extents_ind[1]]*np.sin(input_params['af_params']['azimuthal_location']*np.pi/180),(input_params['af_params']['seperation_distance']+np.sign(input_params['af_params']['seperation_distance'])*D*R)*np.ones(extents_ind[1]-extents_ind[0])]).T,axis = 0)
        af_norms = np.expand_dims(np.array([np.zeros(extents_ind[1]-extents_ind[0]),np.zeros(extents_ind[1]-extents_ind[0]),np.ones(extents_ind[1]-extents_ind[0])]).T,axis = 0)
        
        loads_af = np.array([-dFx_af,np.zeros(dFz_af.shape),-dFz_af]).transpose(1,2,0)


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
    # dFx = np.zeros(dFz.shape)
    loads = np.array([dFy,-dFx,dFz]).transpose(1,2,0)

    # lifting line nodes and normals
    lifting_line_nodes = np.expand_dims(np.array([ac.rotors[0].R*ac.rotors[0].blades[0].r,np.zeros(N_elements),np.zeros(N_elements)]).T,axis = 0)
    lifting_line_norms = np.expand_dims(np.array([np.zeros(N_elements),np.zeros(N_elements),np.ones(N_elements)]).T,axis = 0)

    if input_params['computational_params']['unsteady_loading']:
        saved_params.update({'t':t,'iterations':iterations,'r_elem':ac.rotors[0].blades[0].elems,'r':ac.rotors[0].blades[0].r,"th":ac.rotors[0].blades[0].th,'airfoil':geom_params['airfoil'],'airfoil_points':input_params['computational_params']['airfoil_elements'],'th_tw':th_tw,'TR':TR,'AR':AR,'R':R,'e':e,'c':ac.rotors[0].blades[0].c,'dpsi':dpsi,'dt':dt,'sos':sos,'N_elements':N_elements,'omega':omega,'psi':psi,'x':x,'lam_af':lam_y,'U':U,'s':s,'th0':th0,'CL':CL,'aoa':aoa_eff,'phi_eff':phi_eff,'dCT':dCT,'dCP':dCP,'loads':loads,'lifting_line_nodes':lifting_line_nodes,'lifting_line_norms':lifting_line_norms,'loads_af':loads_af,'af_norms':af_norms,'af_nodes':af_nodes})
    else:
        saved_params.update({'t':t,'iterations':iterations,'r_elem':ac.rotors[0].blades[0].elems,'r':ac.rotors[0].blades[0].r,"th":ac.rotors[0].blades[0].th,'airfoil':geom_params['airfoil'],'airfoil_points':input_params['computational_params']['airfoil_elements'],'th_tw':th_tw,'TR':TR,'AR':AR,'R':R,'e':e,'c':ac.rotors[0].blades[0].c,'dpsi':dpsi,'dt':dt,'sos':sos,'N_elements':N_elements,'omega':omega,'psi':psi,'th0':th0,'CL':ac.rotors[0].blades[0].CL,'CD':ac.rotors[0].blades[0].CD,'aoa':ac.rotors[0].blades[0].aoa,'dCT':dCT,'dCP':dCP,'loads':loads,'lifting_line_nodes':lifting_line_nodes,'lifting_line_norms':lifting_line_norms})

    # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
    # ax.plot(psi,w_af[::25,:,-2].T)
    # ax.plot(psi,aoa_eff[:,-2])
    # fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
    # ax.plot(psi,np.gradient(loads_af[:,-2,-1]))
    # ax.plot(psi,np.gradient(loads_af[:,-2,0]))
    # ax.plot(psi,np.gradient(dFz[:,-20]))
    # # ax.set_xlim([0,np.pi])

