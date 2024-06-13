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


case_name = 'model_360'
case_dir = os.path.join(os.path.dirname(__file__),case_name)
case = caseName(globalFolderName = case_dir,caseNameFile=f'{case_name}.nam')

if os.path.exists(case_dir):
    rmtree(case_dir)
os.mkdir(case_dir)

#%%
N_elements = 48

Nb = 1
R = 38.75/39.37
e = .268*R
N_rotor = 1
c = 3/39.37
th_tw = 0*np.pi/180
th0 = 10*np.pi/180
N_elements = 48
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

dpsi = 1*np.pi/180
iterations = int(2*(2*np.pi)/dpsi)
psi = np.arange(iterations+1)*dpsi

dt = dpsi/omega
t = np.arange(iterations+1)*dt

#%%
geom_dir = os.path.join(os.path.dirname(__file__),'validation','model_360','Boeing360_DegenGeom.csv')
[dataSorted, indHeader] = AnalyzeDegenGeom(geom_dir)
geomParams = ProcessGeom(dataSorted, indHeader, 0.25, Nb, 1)

# nodes = np.array([geomParams['surfNodes'][:,1],geomParams['surfNodes'][:,0],geomParams['surfNodes'][:,-1]]).T
# norms = np.array([geomParams['surfNorms'][:,1],geomParams['surfNorms'][:,0],geomParams['surfNorms'][:,-1]]).T

nodes = geomParams['surfNodes'].reshape(geomParams['pntsPerXsec'],geomParams['nXsecs'],3,order = 'F')
norms = geomParams['surfNorms'].reshape(geomParams['pntsPerXsec'],geomParams['nXsecs'],3,order = 'F')


#%%
atmos = Atmosphere()
ac = aircraft(N_rotor)
ac.rotors = [rotor(Nb = Nb,R = R,e = e,c = c,th0 = th0,th_tw = th_tw,N_elements = N_elements,af = af,Cl_a=Cl_a,origin = origin,omega = omega,V_c = V_c,C_T_target = C_T_target,atmos=atmos) for r_iter in range(ac.N_rotor)]

sol = opt.newton(trim,x0 = ac.rotors[0].th0,args=(ac.rotors[0],ac.rotors[0].blades[0]),tol=5e-6,full_output=False)
ac.rotors[0].th0 = sol
ac.rotors[0].blades[0].set_twist(ac.rotors[0])
ac.rotors[0].blades[0].set_loads()
lam_bemt = ac.rotors[0].blades[0].lam

# dFx = np.ones((iterations+1,Nb,N_elements))*ac.rotors[0].blades[0].dCP*0.5*rho*np.pi*R**2*(omega*R)**3
# dFy = np.zeros((iterations+1,Nb,N_elements))
# dFz = np.ones((iterations+1,Nb,N_elements))*ac.rotors[0].blades[0].dCT*0.5*rho*np.pi*R**2*(omega*R)**2
# dr = np.mean(np.diff(ac.rotors[0].blades[0].r))
dFx= ac.rotors[0].blades[0].dCP*rho*np.pi*R**2*(omega*R)**2
dFz = ac.rotors[0].blades[0].dCT*rho*np.pi*R*(omega*R)**2
dFy  = np.zeros(N_elements)
# dFz  = np.zeros(N_elements)
# dFx  = np.zeros(N_elements)

loads = np.array([dFy,-dFx,dFz]).T
loads = (np.expand_dims(loads,axis = -1)*np.ones(iterations+1)).transpose(-1,0,1)

lifting_line_nodes = np.expand_dims(np.array([ac.rotors[0].R*ac.rotors[0].blades[0].r,np.zeros(N_elements),np.zeros(N_elements)]).T,axis = 0)
lifting_line_norms = np.expand_dims(np.array([np.zeros(N_elements),np.zeros(N_elements),np.ones(N_elements)]).T,axis = 0)

#%%

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

observer_coordinates = np.array([[0,9.238488,9.238488*np.tan(-15*np.pi/180)],[0,9.238488,9.238488*np.tan(-6*np.pi/180)],[0,9.238488,0],[0,9.238488,9.238488*np.tan(6*np.pi/180)]])

# observer_in = ObserverIn(nt = len(psi),Title='mic array',attachedTo = 'aircraft',tMin = t[0],tMax=t[-1],nbx=nbx,xMin=xMin,xMax=xMax,nby=nby,yMin=yMin,yMax=yMax,nbz = nbz,zMin = zMin,zMax=zMax,highPassFrequency = 1)
observer_in = ObserverIn(nt = len(psi),Title='mic array',attachedTo = 'aircraft',tMin = t[0],tMax=t[-1],fileName='observer.ascii',highPassFrequency=0.1)

aircraft_container = ContainerIn(Title='aircraft',nbContainer = 1)
rotor_container = ContainerIn(Title='rotor',nbContainer = 2*Nb,nbBase=1)
rotor_cb = CB(Title='rotation',Rotation = True,AngleType='KnownFunction',Omega=omega,AxisValue=[0,0,1])

nml.extend([environment_in,environment_const,observer_in,aircraft_container,rotor_container,rotor_cb])


for b_iter in range(Nb):
    nml.append(ContainerIn(Title=f'blade {b_iter} loading',nbBase=1,patchGeometryFile = f'lifting_line_geometry.dat',patchLoadingFile = f'loading_blade_{b_iter}.dat'))
    # nml.append(ContainerIn(Title=f'blade {b_iter} loading',nbBase=1,patchGeometryFile = f'lifting_line_geometry.dat',patchLoadingFile = f'constant_loading.dat',dtau = t[-1]/(nt-1)))
    nml.append(CB(Title=f'blade {b_iter} azimuthal offset',AxisValue=[0,0,1],AngleValue=2*np.pi/Nb*b_iter))

    nml.append(ContainerIn(Title=f'blade {b_iter} thickness',nbBase=3,patchGeometryFile = f'blade_geometry.dat',dtau = t[-1]/(nt-1)))
    nml.append(CB(Title=f'blade {b_iter} azimuthal offset',AxisValue=[0,0,1],AngleValue=2*np.pi/Nb*b_iter))
    nml.append(CB(Title=f'blade {b_iter} th0',AxisValue=[1,0,0],AngleValue=ac.rotors[0].th0))
    nml.append(CB(Title=f'Aligns blade span with positive x-direction',AxisValue=[0,0,-1],AngleValue=np.pi/2))

write_nml_file(os.path.join(case_dir,f'{case_name}.nam'),nml)
write_nml_file(os.path.join(os.path.dirname(case_dir),'cases.nam'),[case])
write_observer_file(os.path.join(case_dir,f'observer.ascii'),observer=observer_coordinates)
constant_compact_geometry_write(os.path.join(case_dir,f'blade_geometry.dat'),nodes=nodes,norms=norms,ascii=False)
constant_compact_geometry_write(os.path.join(case_dir,f'lifting_line_geometry.dat'),nodes=lifting_line_nodes,norms=lifting_line_norms,ascii=False)

# constant_compact_loading_write(os.path.join(case_dir,f'constant_loading.dat'), loads =loads,ascii = False)

for b_iter in range(Nb):
    aperiodic_compact_loading_write(os.path.join(case_dir,f'loading_blade_{b_iter}.dat'),t = t, loads = loads,ascii = False)
#     aperiodic_compact_geometry_write(os.path.join(case_dir,f'geometry_blade_{b_iter}.dat'),t = t,nodes=nodes[:,b_iter],norms=norms[:,b_iter],ascii = False)



#%%
# fig,ax = plt.subplots(1,1, figsize = (6.4,4.5),subplot_kw=dict(projection='3d'))
# ax.plot(geomParams2[:,0],geomParams2[:,1],geomParams2[:,2])
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')

# ax.plot_wireframe(blocks[0].X.squeeze(),blocks[0].Y.squeeze(),blocks[0].Z.squeeze())
# ax.scatter(nodes[:,10:12,0],nodes[:,10:12,1],nodes[:,10:12,-1])
