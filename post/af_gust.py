import numpy as np
import matplotlib.pyplot as plt
import aerosandbox as asb

#%%

af = 'naca0012'
c = 0.381/10
M = 0.3
sos = 343
nu = 14.88e-6
th0 = 3*np.pi/180
iterations = 5000
ds = 0.1

Re = M*sos*c/nu

gamma = 0.2*M*sos*c
r_c = 0.05*c
y = -0.26*c
x = (np.arange(iterations+1)/(iterations)*ds/2*iterations-ds/4*iterations)*c
r = np.sqrt(x**2+y**2)

w = gamma/(2*np.pi)*-x/r**2*(1-np.exp(-r**2/r_c**2))
aoa = th0-np.arctan2(w,M*sos)
aoa[-1] = aoa[0]

# Indicial response function coefficients and exponents (these are derived from CFD data and given by Leishman)
A1 = 0.5
b1 = .13
A2 = 0.5
b2 = 1.0

#%%

G1 = 0.67
g1 = 0.1753
G2 = 0.33
g2 = 1.637


U = M*sos
beta = np.sqrt(1-M**2)
s = (x+x[-1])/(c/2)
# ds = np.diff(s[:2])[0]

aoa_eff_1 = np.zeros(iterations+1)
aoa_eff_2 = np.zeros(iterations+1)
X_temp = 0
Y_temp = 0
E1_temp = 0
E2_temp = 0

for i in range(iterations+1):

    X = X_temp*np.exp(-b1*ds)+A1*(aoa[i]-aoa[i-1])*np.exp(-b1*2*ds/2)
    Y = Y_temp*np.exp(-b2*ds)+A2*(aoa[i]-aoa[i-1])*np.exp(-b2*2*ds/2)
    aoa_eff_1[i] = aoa[i]-X-Y
    X_temp = X
    Y_temp = Y

    E1 = E1_temp*np.exp(-g1*(1-M**2)*ds)+G1*(w[i]-w[i-1])*np.exp(-g1*(1-M**2)*2*ds/2)
    E2 = E2_temp*np.exp(-g2*(1-M**2)*ds)+G2*(w[i]-w[i-1])*np.exp(-g2*(1-M**2)*2*ds/2)
    aoa_eff_2[i] = 1/(U*np.sqrt(1-M**2))*(w[i]-E1-E2)
    E1_temp = E1
    E2_temp = E2


af = asb.Airfoil(af)
aero = af.get_aero_from_neuralfoil(alpha=aoa_eff_1*180/np.pi, Re=Re*np.ones(iterations+1),mach =M*np.ones(iterations+1),model_size='xlarge')

daoa = 0.5
aoa_polar = np.arange(20/daoa+1)*daoa-5
aero_polar = af.get_aero_from_neuralfoil(alpha=aoa_polar, Re=Re,mach =M,model_size='xlarge')
a0_ind = np.abs(aero_polar['CL']).argmin()
CL_max_ind = np.abs(np.gradient(aero_polar['CL'])).argmin()
CL_alpha = (aero_polar['CL'][int(0.5*(CL_max_ind+a0_ind))]-aero_polar['CL'][a0_ind])/(aoa_polar[int(0.5*(CL_max_ind+a0_ind))]-aoa_polar[a0_ind])*180/np.pi


fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
ax.plot(x/c+2*ds,aoa*180/np.pi)
ax.plot(x/c+2*ds,aoa_eff_1*180/np.pi)
ax.plot(x/c+2*ds,(-aoa_eff_2+th0)*180/np.pi)
ax.set_xlim([-5,5])

fig,ax = plt.subplots(1,1, figsize = (6.4,4.5))
ax.plot(x/c+2*ds,(aoa_eff_1-aoa_polar[a0_ind])*CL_alpha)
ax.plot(x/c+2*ds,(-aoa_eff_2+th0-aoa_polar[a0_ind])*CL_alpha)
ax.plot(x/c+2*ds,aero['CL'])
ax.set_xlim([-5,5])
# ax.set_ylim([-.25,.1])
ax.grid()
ax.legend(['Linear','Table Lookup'])
