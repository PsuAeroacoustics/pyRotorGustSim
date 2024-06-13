#!/usr/bin/env python3

from geometry import *
from bemt import *
import os
import numpy as np
import matplotlib.pyplot as plt
import aerosandbox as asb
import scipy.optimize as opt
import matplotlib.colors as mcolors
from wopwop_input_generator import *
from shutil import rmtree
from AnalyzeDegenGeom import *
from ProcessGeom import *

cmap = plt.cm.Spectral.reversed()


#%%

case_name = 'single_blade_arbitrary_gust'
case_dir = os.path.join(os.path.dirname(__file__),case_name)

if os.path.exists(case_dir):
    rmtree(case_dir)
os.mkdir(case_dir)

#%% Rotor parameters and operating conditions

N_elements = 48

Nb = 1
R = 38.75/39.37
e = .268*R
N_rotor = 1
c = 3/39.37
th_tw = 0*np.pi/180
th0 = 10*np.pi/180
Cl_a = 2*np.pi
af = asb.Airfoil("naca0009")
origin = [0,0,0]

rho = 1.125
sos = 343
nu = 14.88e-6

# rotational rate [rad/s]
omega = .702*sos/R
V_c = 0
C_T_target = .006/4

dpsi = .25*np.pi/180
iterations = int(2*(2*np.pi)/dpsi)
psi = np.arange(iterations)*dpsi

dt = dpsi/omega
t = np.arange(iterations)*dt


#%% Indicial response function
# # Wagner function approximation
# A1 = 0.165
# b1 = .0455
# A2 = 0.335
# b2 = .3
# Kussner function approximation
# A1 = 0.5
# b1 = .13
# A2 = 0.5
# b2 = 1
# CFD data
A1 = 0.67
b1 = .1753
A2 = 0.33
b2 = 1.637

#%%

atmos = Atmosphere()
ac = aircraft(N_rotor)
ac.rotors = [rotor(Nb = Nb,R = R,e = e,c = c,th0 = th0,th_tw = th_tw,N_elements = N_elements,af = af,Cl_a=Cl_a,origin = origin,omega = omega,V_c = V_c,C_T_target = C_T_target,atmos=atmos) for r_iter in range(ac.N_rotor)]

sol = opt.newton(trim,x0 = ac.rotors[0].th0,args=(ac.rotors[0],ac.rotors[0].blades[0]),tol=5e-6,full_output=False)
ac.rotors[0].th0 = sol
ac.rotors[0].blades[0].set_twist(ac.rotors[0])
ac.rotors[0].blades[0].set_loads()
lam_bemt = ac.rotors[0].blades[0].lam

dFx= ac.rotors[0].blades[0].dCP*rho*np.pi*R**2*(omega*R)**2
dFz = ac.rotors[0].blades[0].dCT*rho*np.pi*R*(omega*R)**2
dFy  = np.zeros(N_elements)
loads = np.array([dFy,-dFx,dFz]).T
loads = (np.expand_dims(loads,axis = -1)*np.ones(iterations+1)).transpose(-1,0,1)


#%%


# gamma_vf = 1.7*np.max(dFz/(ac.rotors[0].blades[0].U*omega*R*rho))/R
gamma_vf = 4.7/R
vf_end_pnts = np.array([[0,0],[0,R]])/R
r_vf = 0.15*c/R
n = 2
gust_ind = int(np.arctan2((vf_end_pnts[-1]-vf_end_pnts[0])[-1],(vf_end_pnts[-1]-vf_end_pnts[0])[0])/dpsi)

h = (np.arange(50+1)*(1.4+.2)/50-.2)/39.37/R



# v_gust = gamma_vf/(2*np.pi*r_vf*(h/r_vf))*(1-np.exp(-1.25643*(h/r_vf)**2))

v_gust = gamma_vf/(2*np.pi*h)*(1-np.exp(-1.25643*(h/r_vf)**2))

v_gust2 = gamma_vf/(2*np.pi)*(h/(r_vf**(2*n)+(h)**(2*n))**(1/n))

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = .15)
ax.plot(h*R*39.37,v_gust*3.281)
ax.plot(h*R*39.37,v_gust2*3.281)
ax.set_ylabel('V [fps]')
ax.set_xlabel('Nozzle Width [in]')
ax.axis([-.2,1.4,0,160])
ax.grid()


# h = np.expand_dims((gust_ind*dpsi-psi)%(2*np.pi),axis = -1)*ac.rotors[0].blades[0].r
# h = np.expand_dims((gust_ind*dpsi-psi)%(2*np.pi),axis = -1)*ac.rotors[0].blades[0].r
h = np.expand_dims(((psi%(2*np.pi)-psi[gust_ind])%(2*np.pi)),axis = -1)*ac.rotors[0].blades[0].r

# v_gust = gamma_vf/(2*np.pi*r_vf*(h/r_vf))*(1-np.exp(-1.25643*(h/r_vf)**2))
# v_gust = gamma_vf/(2*np.pi*h)*(1-np.exp(-1.25643*(h/r_vf)**2))
v_gust = gamma_vf/(2*np.pi)*(h/(r_vf**(2*n)+(h)**(2*n))**(1/n))

nan_ind = np.where(np.isnan(v_gust))
v_gust[nan_ind] = v_gust[(nan_ind[0]+1,nan_ind[1])]
lam_gust = v_gust/(omega*R)

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = .15)
ax.plot(psi*180/np.pi,v_gust[:,-20]*3.281)
ax.grid()

# fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
# plt.subplots_adjust(left = .15)
# ax.plot(h[:,-1],v_gust[:,-1])
# ax.grid()
# h = ((np.arange(50)*1.4/50)/39.37)/R
# # v_ind_2 = 1/(2*np.pi)*gamma_vf*h/(r_vf**(2*n)+h**(2*n))**(1/n)

# h = np.linalg.norm(((vf_end_pnts[-1]-vf_end_pnts[0])-np.array((np.cos(psi[:360]),np.sin(psi[:360]))).T),axis = -1)

# v_ind = gamma_vf/(2*np.pi*r_vf*(h/r_vf))*(1-np.exp(-1.25643*(h/r_vf)**2))


# fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
# plt.subplots_adjust(left = .15)
# ax.plot(psi[:360]*180/np.pi,v_ind)
# ax.set_ylabel('$V$')
# ax.set_xlabel('Nozzle width')
# # ax.set_xlim([-0.2,1.4])
# # ax.set_ylim([0,160])
# ax.grid()


#%%


# lam = np.ones((iterations,N_elements))*lam_bemt
lam = lam_bemt+lam_gust
U = np.sqrt(lam**2+ac.rotors[0].blades[0].r**2)
beta = np.sqrt(1-(U*omega*R/sos)**2)
phi = np.arctan2(lam,ac.rotors[0].blades[0].r)
aoa = ac.rotors[0].blades[0].th-phi

s = omega*R*ac.rotors[0].blades[0].r*np.expand_dims(t,axis = -1)/(c/2)
ds = np.diff(s,axis = 0)[0]

aoa_eff = np.zeros((iterations,N_elements))
X_temp = np.zeros((N_elements))
Y_temp = np.zeros((N_elements))

for i in range(iterations):

    X = X_temp*np.exp(-b1*beta[i]**2*ds)+A1*omega*R*(lam_gust[i]-lam_gust[i-1])*np.exp(-b1*beta[i]**2*ds)**(1/2)
    Y = Y_temp*np.exp(-b2*beta[i]**2*ds)+A2*omega*R*(lam_gust[i]-lam_gust[i-1])*np.exp(-b2*beta[i]**2*ds)**(1/2)
    dCL = 2*np.pi/(beta[i]*U[i]*omega*R)*(lam_gust[i]*omega*R-X-Y)
    aoa_eff[i] = dCL/(2*np.pi)

    # X = X_temp*np.exp(-b1*ds)+A1*(aoa[i]-aoa[i-1])*np.exp(-b1*ds/2)
    # Y = Y_temp*np.exp(-b2*ds)+A2*(aoa[i]-aoa[i-1])*np.exp(-b2*ds/2)
    # aoa_eff[i] = aoa[i]-X-Y

    X_temp = X
    Y_temp = Y

# aoa_eff = aoa_eff+ac.rotors[0].blades[0].aoa
phi_eff = ac.rotors[0].blades[0].th - aoa_eff

Re = np.ones(aoa_eff.shape)*ac.rotors[0].blades[0].Re
M = np.ones(aoa_eff.shape)*ac.rotors[0].blades[0].M

CL,CD = get_af_coeffs(ac.rotors[0].blades[0].af,aoa_eff*180/np.pi,Re,M)

# CL = 2*np.pi*aoa_eff+ac.rotors[0].blades[0].aoa
# CD = np.zeros(CL.shape)

dCz = CL*np.cos(phi_eff)-CD*np.sin(phi_eff)
dCx = CL*np.sin(phi_eff)+CD*np.cos(phi_eff)

dCT = 0.5*ac.rotors[0].c/(np.pi*ac.rotors[0].R)*U**2*dCz
dCP = 0.5*ac.rotors[0].c/(np.pi*ac.rotors[0].R)*ac.rotors[0].blades[0].r*U**2*dCx

dFx= dCP*atmos.rho*np.pi*ac.rotors[0].R**2*(ac.rotors[0].omega*ac.rotors[0].R)**2
dFz = dCT*atmos.rho*np.pi*ac.rotors[0].R*(ac.rotors[0].omega*ac.rotors[0].R)**2

# dFz = 0.5*atmos.rho*(U*omega*R)**2*ac.rotors[0].c*dCz
# dFx = 0.5*atmos.rho*(U*omega*R)**2*ac.rotors[0].c*dCx*ac.rotors[0].blades[0].r*ac.rotors[0].R*ac.rotors[0].omega
dFy = np.zeros(dFz.shape)
# dFx = np.zeros(dFz.shape)
# loads = np.array([dFy[0],-dFx[0],dFz[0]])

loads = np.array([dFy,-dFx,dFz]).transpose(1,2,0)



#%% Process blade geometry

geom_dir = os.path.join(os.path.dirname(__file__),'validation','model_360','Boeing360_DegenGeom.csv')
[dataSorted, indHeader] = AnalyzeDegenGeom(geom_dir)
geomParams = ProcessGeom(dataSorted, indHeader, 0.25, Nb, 1)

nodes = geomParams['surfNodes'].reshape(geomParams['pntsPerXsec'],geomParams['nXsecs'],3,order = 'F')
norms = geomParams['surfNorms'].reshape(geomParams['pntsPerXsec'],geomParams['nXsecs'],3,order = 'F')

lifting_line_nodes = np.expand_dims(np.array([ac.rotors[0].R*ac.rotors[0].blades[0].r,np.zeros(N_elements),np.zeros(N_elements)]).T,axis = 0)
lifting_line_norms = np.expand_dims(np.array([np.zeros(N_elements),np.zeros(N_elements),np.ones(N_elements)]).T,axis = 0)


#%%

nml = []

environment_in = EnvironmentIn(debugLevel = 4,ASCIIOutputFlag=True,totalNoiseFlag=True,pressureFolderName = '/')
environment_const = EnvironmentConstants()

nt = int(iterations/2)


# radius = 3*R
# nbTheta = 18
# nbPsi = 1
# thetaMin = 90*np.pi/180
# thetaMax = 270*np.pi/180
# psiMin = 50.2*np.pi/180
# psiMax = 50.2*np.pi/180

radius = 3*R
# th = np.arange(15)*(2*np.pi)/15
# th = np.array([104-30+90,107-30+90])*np.pi/180
th = (np.array([66,76,106,118,132])+60)*np.pi/180

# th = (np.array([61.6,75.5,91.4,107,120.9,133,148])+60)*np.pi/180
phi = 50.2*np.pi/180

observer_coordinates = np.array([radius*np.cos(th),radius*np.sin(th),radius*np.sin(phi)*np.ones(len(th))]).T
observer_in = ObserverIn(nt = len(psi),Title='mic array',attachedTo = 'aircraft',tMin = t[int(iterations/2)],tMax=t[-1],fileName='observer.ascii',highPassFrequency=0.1)
write_observer_file(os.path.join(case_dir,f'observer.ascii'),observer=observer_coordinates)

# observer_coordinates = np.array([[-9.238488,0,9.238488*np.tan(-15*np.pi/180)],[-9.238488,0,9.238488*np.tan(-6*np.pi/180)],[-9.238488,0,0],[-9.238488,0,9.238488*np.tan(6*np.pi/180)]])
# observer_in = ObserverIn(nt = len(psi),Title='mic array',attachedTo = 'aircraft',tMin = t[0],tMax=t[-1],nbx=nbx,xMin=xMin,xMax=xMax,nby=nby,yMin=yMin,yMax=yMax,nbz = nbz,zMin = zMin,zMax=zMax,highPassFrequency = 1)
# observer_in = ObserverIn(nt = len(psi),Title='mic array',attachedTo = 'aircraft',tMin = t[0],tMax=t[-1],fileName='observer.ascii',highPassFrequency=0.1)
# write_observer_file(os.path.join(case_dir,f'observer.ascii'),observer=observer_coordinates)

# observer_in = ObserverIn(nt = len(psi),Title='mic array',attachedTo = 'aircraft',radius=radius,nbTheta=nbTheta,nbPsi=nbPsi,thetaMin=thetaMin,thetaMax=thetaMax,psiMin=psiMin,psiMax=psiMax,tMin = t[int(iterations/2)],tMax=t[-1],highPassFrequency=0.1)

aircraft_container = ContainerIn(Title='aircraft',nbContainer = 1)
rotor_container = ContainerIn(Title='rotor',nbContainer = 2*Nb,nbBase=1)
rotor_cb = CB(Title='rotation',Rotation = True,AngleType='KnownFunction',Omega=omega,AxisValue=[0,0,1])
nml.extend([environment_in,environment_const,observer_in,aircraft_container,rotor_container,rotor_cb])


for b_iter in range(Nb):
    nml.append(ContainerIn(Title=f'blade {b_iter} loading',nbBase=1,patchGeometryFile = f'lifting_line_geometry.dat',patchLoadingFile = f'loading_blade_{b_iter}.dat',dtau = t[-1]/(nt-1)))
    # nml.append(ContainerIn(Title=f'blade {b_iter} loading',nbBase=1,patchGeometryFile = f'lifting_line_geometry.dat',patchLoadingFile = f'constant_loading.dat',dtau = t[-1]/(nt-1)))
    nml.append(CB(Title=f'blade {b_iter} azimuthal offset',AxisValue=[0,0,1],AngleValue=2*np.pi/Nb*b_iter))

    nml.append(ContainerIn(Title=f'blade {b_iter} thickness',nbBase=3,patchGeometryFile = f'blade_geometry.dat',dtau = t[-1]/(nt-1)))
    nml.append(CB(Title=f'blade {b_iter} azimuthal offset',AxisValue=[0,0,1],AngleValue=2*np.pi/Nb*b_iter))
    nml.append(CB(Title=f'blade {b_iter} th0',AxisValue=[1,0,0],AngleValue=ac.rotors[0].th0))
    nml.append(CB(Title=f'Aligns blade span with positive x-direction',AxisValue=[0,0,-1],AngleValue=np.pi/2))

write_nml_file(os.path.join(case_dir,f'{case_name}.nam'),nml)

case = caseName(globalFolderName = case_dir,caseNameFile=f'{case_name}.nam')
write_nml_file(os.path.join(os.path.dirname(case_dir),'cases.nam'),[case])


constant_compact_geometry_write(os.path.join(case_dir,f'blade_geometry.dat'),nodes=nodes,norms=norms,ascii=False)
constant_compact_geometry_write(os.path.join(case_dir,f'lifting_line_geometry.dat'),nodes=lifting_line_nodes,norms=lifting_line_norms,ascii=False)
# constant_compact_loading_write(os.path.join(case_dir,f'constant_loading.dat'), loads =loads,ascii = False)

# periodic_compact_loading_write(os.path.join(case_dir,f'loading_blade_{b_iter}.dat'),keys = 180/np.pi*psi[:360] ,period = 2*np.pi/omega, loads=loads[:360],ascii = False)

for b_iter in range(Nb):
    aperiodic_compact_loading_write(os.path.join(case_dir,f'loading_blade_{b_iter}.dat'),t = t, loads = loads,ascii = False)


#%%

# fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
# plt.subplots_adjust(left = .2)
# ax.plot(psi*180/np.pi,np.gradient(aoa[:,30]))
# ax.plot(psi*180/np.pi,np.gradient(aoa_eff[:,30]))
# ax.set_ylabel('$\partial \\alpha / \partial \psi$')
# ax.set_xlabel('$\psi \ [deg]$')
# # ax.set_xlim([0,360])
# # ax.set_ylim([-0.01,0.01])
# ax.legend(['w/o indicial response','w/ indicial response'])
# ax.grid()

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
ax.plot(psi*180/np.pi,aoa[:,30]*180/np.pi)
ax.plot(psi*180/np.pi,aoa_eff[:,30]*180/np.pi)
ax.set_ylabel('$ \\alpha \ [deg]$')
ax.set_xlabel('$\psi \ [deg]$')
# ax.set_xlim([0,360])
# ax.set_ylim([-0.01,0.01])
ax.legend(['w/o indicial response','w/ indicial response'])
ax.grid()


fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = .15)
ax.plot(psi/(2*np.pi),dCT[:,30])
ax.set_ylabel('$\partial CT/\partial \psi$')
ax.set_xlabel('Rotation')
ax.grid()


# fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
# ax.plot(ac.rotors[0].blades[0].r,ac.rotors[0].blades[0].lam)
# ax.set_xlabel('r/R')
# ax.grid()
# ax.set_ylabel('$\lambda$')
# ax.set_xlim([0,1])
# ax.set_ylim([0,.1])


# fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
# ax.plot(ac.rotors[0].blades[0].r,ac.rotors[0].blades[0].dCT)
# ax.set_xlabel('r/R')
# ax.grid()
# # ax.set_ylabel('$\dC_T$')
# ax.set_xlim([0,1])
# ax.set_ylim([0,.03])

# dcm = np.array([[np.cos(psi),np.sin(psi),np.zeros(len(psi))],[-np.sin(psi),np.cos(psi),np.zeros(len(psi))],[np.zeros(len(psi)),np.zeros(len(psi)),np.ones(len(psi))]])

# V_ind_trans = np.zeros(V_ind.shape).transpose(2,0,1,-1)
# for b_iter in range(Nb):
#     dcm(psi+(2*np.pi/Nb*b_iter))
#     V_ind_trans[b_iter] = np.matmul(dcm(psi+(2*np.pi/Nb*b_iter)),V_ind[:,:,b_iter])


# V_th = (-V_ind[:360,0,0].T*np.sin(psi[:360])+V_ind[:360,1,0].T*np.cos(psi[:360])).T
# V_r = (V_ind[:360,1,0].T*np.sin(psi[:360])+V_ind[:360,0,0].T*np.cos(psi[:360])).T

# np.matmul(dcm,V_ind[:,:,0])[:10,1,30]

# levels = np.linspace(np.min(V_ind[:90,-1,0]), np.max(V_ind[:90,-1,0]), 50)

# fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
# quant = lam[:360]
# levels = np.linspace(np.min(quant), np.max(quant), 50)
# dist = ax.contourf(psi[:360], ac.rotors[0].blades[0].r, quant.T, levels=levels,cmap = cmap,norm=mcolors.CenteredNorm())
# ax.set_ylim(0, 1)
# ax.set_yticks(ax.get_yticks()[::2])
# cbar = fig.colorbar(dist,format = '%1.2e',pad = .075)
# cbar.ax.set_ylabel('$\lambda$')

fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
quant = dCT[:int(360/dpsi)]
levels = np.linspace(np.min(quant), np.max(quant), 50)
dist = ax.contourf(psi[:int(360/dpsi)], ac.rotors[0].blades[0].r, quant.T, levels=levels,cmap = cmap,norm=mcolors.CenteredNorm())
# ax.set_ylim(0, 1)
ax.set_yticks(ax.get_yticks()[::2])
cbar = fig.colorbar(dist,format = '%1.2e',pad = .075)
cbar.ax.set_ylabel('$\partial CT/\partial \psi$')

fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
quant = np.gradient(dCT[:int(360/dpsi)],axis = 0)
levels = np.linspace(np.min(quant), np.max(quant), 50)
dist = ax.contourf(psi[:int(360/dpsi)], ac.rotors[0].blades[0].r, quant.T, levels=levels,cmap = cmap,norm=mcolors.CenteredNorm())
ax.set_ylim(0, 1)
ax.set_yticks(ax.get_yticks()[::2])
cbar = fig.colorbar(dist,format = '%1.2e',pad = .075)
cbar.ax.set_ylabel('$\partial CT/\partial \psi$')

# fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
# quant = lam[:360]
# levels = np.linspace(np.min(quant), np.max(quant), 50)
# dist = ax.contourf(psi[:360], ac.rotors[0].blades[0].r, quant.T, levels=levels,cmap = cmap,norm=mcolors.CenteredNorm())
# ax.set_ylim(0, 1)
# ax.set_yticks(ax.get_yticks()[::2])
# cbar = fig.colorbar(dist,format = '%1.2e',pad = .075)
# cbar.ax.set_ylabel('$\lambda_i$')



