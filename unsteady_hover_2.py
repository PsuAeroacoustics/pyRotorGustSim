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

    #     return err

def compute_induced_velocity(vf_coord,b_coord,gamma_vf,r_vf):

    b_coord_shape = np.array(b_coord.shape)
    b_coord = b_coord.reshape((3,np.prod(b_coord_shape[np.invert(b_coord_shape==3)])))

    r1 = (b_coord.T-vf_coord[0]).T
    r2 = (b_coord.T-vf_coord[-1]).T
    r1_norm = r1/np.linalg.norm(r1,axis = 0)
    r2_norm = r2/np.linalg.norm(r2,axis = 0)
    l_v = vf_coord[-1]-vf_coord[0]
    l_v_norm = l_v/np.linalg.norm(l_v)

    n_sin = np.cross(l_v_norm,r1_norm,axis = 0)
    e = n_sin/np.linalg.norm(n_sin,axis = 0)
    h = np.linalg.norm(r1*(n**-1*n_sin),axis = 0)
    V_ind = (1/(4*np.pi)*gamma_vf*h/(r_vf**(2*n)+h**(2*n))**(1/n)*(np.dot(l_v_norm,r1_norm-r2_norm)*e)).reshape(b_coord_shape)

    return V_ind

#%%

case_name = 'stationary_vortex_.5R'
case_dir = os.path.join(os.path.dirname(__file__),case_name)

if os.path.exists(case_dir):
    rmtree(case_dir)
os.mkdir(case_dir)

#%% Rotor parameters and operating conditions

N_elements = 48

Nb = 4
R = 1.539748
e = .268*R
N_rotor = 1
c = 0.134239
th_tw = -9.3*np.pi/180
th0 = 10*np.pi/180
N_elements = 48
Cl_a = 2*np.pi
af = asb.Airfoil("vr12")
origin = [0,0,0]

rho = 1.125
sos = 343
nu = 14.88e-6

# rotational rate [rad/s]
omega = .636*sos/R
V_c = 0
C_T_target = .0701*Nb*c/(np.pi*R)

dpsi = 1*np.pi/180
iterations = int(2*(2*np.pi)/dpsi)
psi = np.arange(iterations+1)*dpsi

dt = dpsi/omega
t = np.arange(iterations+1)*dt

#%% Vortex parameters

vf_end_pnts = np.array([[0,0,1.01*R],[0,R,1.01*R]])/R
V_vf = np.array([0,0,0*R/(omega/(2*np.pi))**-1/2])/R
r_vf = 0.05*c/R
n = 2

#%% Indicial response function
# # Wagner function approximation
# A1 = 0.165
# b1 = .0455
# A2 = 0.335
# b2 = .3
# Kussner function approximation
A1 = 0.5
b1 = .13
A2 = 0.5
b2 = 1

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
 # sets vortex strength to be the maximum of bound circulation
gamma_vf = -np.max(dFz/(ac.rotors[0].blades[0].U*omega*R*rho))/R

#%% Process blade geometry

geom_dir = os.path.join(os.path.dirname(__file__),'validation','model_360','Boeing360_DegenGeom.csv')
[dataSorted, indHeader] = AnalyzeDegenGeom(geom_dir)
geomParams = ProcessGeom(dataSorted, indHeader, 0.25, Nb, 1)

nodes = geomParams['surfNodes'].reshape(geomParams['pntsPerXsec'],geomParams['nXsecs'],3,order = 'F')
norms = geomParams['surfNorms'].reshape(geomParams['pntsPerXsec'],geomParams['nXsecs'],3,order = 'F')

lifting_line_nodes = np.expand_dims(np.array([ac.rotors[0].R*ac.rotors[0].blades[0].r,np.zeros(N_elements),np.zeros(N_elements)]).T,axis = 0)
lifting_line_norms = np.expand_dims(np.array([np.zeros(N_elements),np.zeros(N_elements),np.ones(N_elements)]).T,axis = 0)


#%%

V_ind = np.zeros((iterations+1,3,ac.rotors[0].Nb,ac.rotors[0].N_elements))
dCT = np.zeros((iterations+1,ac.rotors[0].Nb,ac.rotors[0].N_elements))

dcm = lambda psi: np.array([[np.cos(psi),np.sin(psi),np.zeros(len(psi))],[-np.sin(psi),np.cos(psi),np.zeros(len(psi))],[np.zeros(len(psi)),np.zeros(len(psi)),np.ones(len(psi))]]).transpose(-1,0,1).squeeze()

# np.matmul(dcm(psi[i]+np.arange(Nb)*(2*np.pi/Nb)),np.array([ac.rotors[0].blades[0].r,np.zeros(N_elements),np.zeros(N_elements)]))
for i in range(iterations+1):

    b_psi = psi[i]+np.arange(Nb)*(2*np.pi/Nb)
    b_coord = np.array((ac.rotors[0].blades[0].r*np.expand_dims(np.cos(b_psi),axis = -1),
                        ac.rotors[0].blades[0].r*np.expand_dims(np.sin(b_psi),axis = -1),
                        np.zeros((Nb,N_elements))))
    
    # b_coord = np.matmul(dcm(psi[i]+np.arange(Nb)*(2*np.pi/Nb)).T,np.array([ac.rotors[0].blades[0].r,np.zeros(N_elements),np.zeros(N_elements)]))

    V_ind_temp = compute_induced_velocity(vf_end_pnts,b_coord,gamma_vf,r_vf)
    # resolves velocity vectors into the rotating frame
    V_ind[i] =np.matmul(dcm(b_psi),V_ind_temp.transpose(1,0,-1)).transpose(1,0,-1)
    vf_end_pnts = vf_end_pnts+V_vf*dt


lam_vf = V_ind[:,-1]/(ac.rotors[0].omega*ac.rotors[0].R)
lam_p = lam_vf+lam_bemt
lam_t = V_ind[:,1]/(ac.rotors[0].omega*ac.rotors[0].R)

U = np.sqrt(lam_p**2+ac.rotors[0].blades[0].r**2+lam_t**2)
phi = np.arctan2(lam_p,lam_t+ac.rotors[0].blades[0].r)
aoa = ac.rotors[0].blades[0].th-phi

s = omega*R*ac.rotors[0].blades[0].r*np.expand_dims(t,axis = -1)/(c/2)

aoa_eff = np.zeros(aoa.shape)
X_temp = np.zeros((Nb,N_elements))
Y_temp = np.zeros((Nb,N_elements))

for i in range(iterations):
    X = X_temp*np.exp(-b1*(s[i+1]-s[i]))+A1*(aoa[i+1]-aoa[i])*np.exp(-b1*(s[i+1]-s[i])/2)
    Y = Y_temp*np.exp(-b2*(s[i+1]-s[i]))+A2*(aoa[i+1]-aoa[i])*np.exp(-b2*(s[i+1]-s[i])/2)
    aoa_eff[i] = aoa[i]-X-Y
    X_temp = X
    Y_temp = Y

Re = np.ones(aoa.shape)*ac.rotors[0].blades[0].Re
M = np.ones(aoa.shape)*ac.rotors[0].blades[0].M

CL,CD = get_af_coeffs(ac.rotors[0].blades[0].af,aoa*180/np.pi,Re,M)

dCz = CL*np.cos(phi)-CD*np.sin(phi)
dCx = CL*np.sin(phi)+CD*np.cos(phi)

dCT = 0.5*ac.rotors[0].c/(np.pi*ac.rotors[0].R)*U**2*dCz
dCP = 0.5*ac.rotors[0].c/(np.pi*ac.rotors[0].R)*ac.rotors[0].blades[0].r*U**2*dCx

dFx= dCP*atmos.rho*np.pi*ac.rotors[0].R**2*(ac.rotors[0].omega*ac.rotors[0].R)**2
dFz = dCT*atmos.rho*np.pi*ac.rotors[0].R*(ac.rotors[0].omega*ac.rotors[0].R)**2

# dFz = 0.5*atmos.rho*(U*omega*R)**2*ac.rotors[0].c*dCz
# dFx = 0.5*atmos.rho*(U*omega*R)**2*ac.rotors[0].c*dCx*ac.rotors[0].blades[0].r*ac.rotors[0].R*ac.rotors[0].omega
dFy = np.zeros(dFz.shape)

# T = np.sum(np.trapz(dFz,x = ac.rotors[0].blades[0].r*R,axis = -1),axis = -1)
# P = np.sum(np.trapz(dFx,x = ac.rotors[0].blades[0].r*R,axis = -1),axis = -1)

# dCT = ac.rotors[0].sigma/(2*ac.rotors[0].Nb)*U**2*dCz
# CT = np.sum(np.trapz(dCT,x = ac.rotors[0].blades[0].r,axis = -1),axis = -1)

# dCP = ac.rotors[0].sigma/(2*ac.rotors[0].Nb)*ac.rotors[0].blades[0].r*U**2*dCx
# CP = np.sum(np.trapz(dCP,x = ac.rotors[0].blades[0].r,axis = -1),axis = -1)

# combines loads into a single matrix of size (Nb x iterations+1 x N_elements x 3)
loads = np.array([dFy,-dFx,dFz]).transpose(2,1,-1,0)

#%% Writes all input files for wopwop


nml = []

environment_in = EnvironmentIn(debugLevel = 4,ASCIIOutputFlag=True,totalNoiseFlag=True,pressureFolderName = '/')
environment_const = EnvironmentConstants()

nt = len(psi)

nbx = 1
xMin = 0
xMax = 0
nby = 1
yMin = 9.238488
yMax = 9.238488
nbz = 5
zMin = -2.475445398884618
zMax = 2.475445398884618

observer_coordinates = np.array([[-9.238488,0,9.238488*np.tan(-15*np.pi/180)],[-9.238488,0,9.238488*np.tan(-6*np.pi/180)],[-9.238488,0,0],[-9.238488,0,9.238488*np.tan(6*np.pi/180)]])
# observer_in = ObserverIn(nt = len(psi),Title='mic array',attachedTo = 'aircraft',tMin = t[0],tMax=t[-1],nbx=nbx,xMin=xMin,xMax=xMax,nby=nby,yMin=yMin,yMax=yMax,nbz = nbz,zMin = zMin,zMax=zMax,highPassFrequency = 1)
observer_in = ObserverIn(nt = len(psi),Title='mic array',attachedTo = 'aircraft',tMin = t[0],tMax=t[-1],fileName='observer.ascii',highPassFrequency=0.1)
write_observer_file(os.path.join(case_dir,f'observer.ascii'),observer=observer_coordinates)

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

for b_iter in range(Nb):
    aperiodic_compact_loading_write(os.path.join(case_dir,f'loading_blade_{b_iter}.dat'),t = t, loads = loads[b_iter],ascii = False)


#%%

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = .15)
ax.plot(psi/(2*np.pi),np.gradient(aoa[:,0,30],axis =0))
ax.plot(psi/(2*np.pi),np.gradient(aoa_eff[:,0,30],axis =0))
ax.set_ylabel('$\partial CT/\partial \psi$')
ax.set_xlabel('Rotation')
ax.grid()

# fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
# plt.subplots_adjust(left = .15)
# ax.plot(psi/(2*np.pi),dCT[:,0,30])
# ax.set_ylabel('$\partial CT/\partial \psi$')
# ax.set_xlabel('Rotation')
# ax.grid()


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
# quant = np.gradient(dCT[:360,0],axis = 0)
# levels = np.linspace(np.min(quant), np.max(quant), 50)
# dist = ax.contourf(psi[:360], ac.rotors[0].blades[0].r, quant.T, levels=levels,cmap = cmap,norm=mcolors.CenteredNorm())
# ax.set_ylim(0, 1)
# ax.set_yticks(ax.get_yticks()[::2])
# cbar = fig.colorbar(dist,format = '%1.2e',pad = .075)
# cbar.ax.set_ylabel('$\partial CT/\partial \psi$')

# fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
# quant = V_ind[:360,-1,0]/(omega*R)
# levels = np.linspace(np.min(quant), np.max(quant), 50)
# dist = ax.contourf(psi[:360], ac.rotors[0].blades[0].r, quant.T, levels=levels,cmap = cmap,norm=mcolors.CenteredNorm())
# ax.set_ylim(0, 1)
# ax.set_yticks(ax.get_yticks()[::2])
# cbar = fig.colorbar(dist,format = '%1.2e',pad = .075)
# cbar.ax.set_ylabel('$\lambda_i$')



