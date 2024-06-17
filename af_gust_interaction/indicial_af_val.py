import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
#%%

cfd_data_f = 'cfd_data.csv'
cfd_data = np.loadtxt(cfd_data_f,delimiter =',').astype(float)

#%%
c = 1
sos = 340
M = .65
gamma = 0.2*sos*M*c
y_v = -.26*c
dist = 40*c
r_c = 0.05*c
n = 1

#%%

dt = c/(8*M*sos)
N_steps = int((dist/(M*sos))/dt)

t = np.arange(N_steps)*dt

x_b,y_b= .25*c,0

x_v =M*sos*t-dist/2
y_v = np.ones(N_steps)*y_v

# r = np.array([x_b,y_b])-np.array([x_v,y_v]).T
# d = np.linalg.norm(r,axis = -1)

# dx = x_b-x_v
# dy = y_b-y_v

# np.linalg.norm((dx,dy),axis = 0)*np.exp(np.arctan2(dy,dx))

# v = 1/(2*np.pi)*gamma/r_c*(d/r_c)/(r_c**(2*n)+d**(2*n))**(1/n)

# d =np.sign(x_b-x_v)*np.sqrt((x_b-x_v)**2+(y_b-y_v)**2)
# v = -1/(2*np.pi)*gamma*d/(r_c**(2*n)+d**(2*n))**(1/n)

# v = gamma/(2*np.pi)*(d/(r_c**(2*n)+d**(2*n))**(1/n))
d = np.sqrt((x_b-x_v)**2+(y_b-y_v)**2)
v = -gamma/(2*np.pi)*(x_b-x_v)/d**2*(1-np.exp(-d**2/r_c**2))

beta = np.sqrt(1-M**2)
s = M*sos*t/(c/2)
ds = s[1]-s[0]

#%%

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = .175)
ax.plot(x_v[1:]/c,(v/(M*sos))[1:])
ax.set_xlabel('$x_v/c$')
ax.grid()
ax.set_ylabel('$aoa$')


#%%
# CFD data
A1 = 0.67
b1 = .1753
A2 = 0.33
b2 = 1.637

dCL = np.zeros(N_steps)
dCL2 = np.zeros(N_steps)

# X_temp,Y_temp =-0.32944852495635807,-0.037184955063786156
X_temp,Y_temp =0,0

for i in range(1,N_steps):
    X = X_temp*np.exp(-b1*beta**2*ds)+A1*(v[i]-v[i-1])*np.exp(-b1*beta**2*ds)**(1/2)
    Y = Y_temp*np.exp(-b2*beta**2*ds)+A2*(v[i]-v[i-1])*np.exp(-b2*beta**2*ds)**(1/2)
    dCL[i] = 2*np.pi/(beta*M*sos)*(v[i]-X-Y)
    X_temp = X
    Y_temp = Y
    if i== 120:
        print(X)
        print(Y)

A1 = 0.527
b1 = .1
A2 = 0.473
b2 = 1.367
X_temp,Y_temp =-0.22907751174485896,-0.06209302379075392

X_temp,Y_temp =0,0
for i in range(1,N_steps):
    X = X_temp*np.exp(-b1*beta**2*ds)+A1*(v[i]-v[i-1])*np.exp(-b1*beta**2*ds)**(1/2)
    Y = Y_temp*np.exp(-b2*beta**2*ds)+A2*(v[i]-v[i-1])*np.exp(-b2*beta**2*ds)**(1/2)
    dCL2[i] = 2*np.pi/(beta*M*sos)*(v[i]-X-Y)
    X_temp = X
    Y_temp = Y
    if i==120:
        print(X)
        print(Y)

#%%

#%%

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = .175)
ax.scatter(cfd_data[:,0],cfd_data[:,-1],c = 'black')
ax.plot(x_v[1:]/c,dCL[1:])
ax.plot(x_v[1:]/c,dCL2[1:],linestyle = '--')
ax.set_xlabel('$x_v/c$')
ax.grid()
ax.set_ylabel('$CL$')
ax.legend(['CFD','Indicial (Linear)', 'Indicial (CFD)'])
ax.axis([-5,5,-.25,.1])