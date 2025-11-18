import numpy as np 
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ["Times New Roman"]
plt.rcParams['font.size'] = 12

#%%

def get_gust_position(**kwargs):
    y_g_func = lambda x: np.tan(kwargs['th_gust'])*(x-np.cos(kwargs['psi_gust_intersect']))+np.sin(kwargs['psi_gust_intersect'])
    gust_ind = np.abs(y_g_func(x_b[:int(psi_max/dpsi)])-y_b[:int(psi_max/dpsi)]).argmin(axis = 0)
    x_g = r_elem*np.cos(psi[gust_ind])
    y_g = r_elem*np.sin(psi[gust_ind])
    psi_g = np.arctan2(y_g,x_g)
    gamma = np.abs(psi_g-kwargs['th_gust'])
    M_tr = kwargs['M_b']/np.sin(gamma)
    if kwargs['th_gust']<=kwargs['psi_gust_intersect']:
        th_rad = np.pi-(np.pi/2-np.arcsin(1/M_tr[M_tr>1])-kwargs['th_gust'])
    else:
        th_rad = (np.pi/2-np.arcsin(1/M_tr[M_tr>1]))+kwargs['th_gust']
    return {'x_g':x_g,'y_g':y_g,'gamma':gamma,'M_tr':M_tr,'th_rad':th_rad}

#%%

R = 0.381
r_c = 0.268
N_elements = 48
sos = 340
dpsi = 1*np.pi/180
omega = np.expand_dims(sos/R*(np.arange(5)*(.8-.4)/4+.4),axis = -1)

r = np.arange(N_elements+1)*(1-r_c)/N_elements+r_c
r_elem = 0.5*(r[1:]+r[:-1])
psi = np.arange(int(2*np.pi/dpsi)+1)*dpsi
x_b = r_elem*np.cos(np.expand_dims(psi,axis = -1))
y_b = r_elem*np.sin(np.expand_dims(psi,axis = -1))
M_b = omega*r_elem*R/sos

#%%

th_gust = 45/2*np.pi/180
psi_gust_intersect = 45*np.pi/180
psi_max = np.pi/2

gust_params = list(map(lambda x: get_gust_position(**{'th_gust':th_gust,'psi_gust_intersect':psi_gust_intersect,'M_b':x}),M_b))

#%%

fig,ax = plt.subplots(1,1, figsize = (4.5,4.5))
plt.subplots_adjust(left = 0.175,right = .95)
ax.plot(np.cos(psi),np.sin(psi))
ax.plot(gust_params[0]['x_g'],gust_params[0]['y_g'])
ax.set_ylabel('$y/R$')
ax.set_xlabel('$x/R$')
ax.grid()

#%%
leg_lab = [f'$M_T = {M}$'for M in np.round(omega.squeeze()*R/sos,2)]

fig,ax = plt.subplots(1,1, figsize = (4.5,3.375))
plt.subplots_adjust(bottom = 0.15,left= 0.15)
for i in range(len(gust_params)):
    ax.scatter(r_elem[gust_params[i]['M_tr']>1],gust_params[i]['M_tr'][gust_params[i]['M_tr']>1])
ax.set_ylabel('$M_{tr}$')
ax.set_xlabel('$r/R$')
ax.legend(leg_lab)
ax.grid()

#%%
fig,ax = plt.subplots(1,1, figsize = (4.5,3.375))
plt.subplots_adjust(bottom = 0.15,left= 0.15)
for i in range(len(gust_params)):
    ax.scatter(r_elem[gust_params[i]['M_tr']>1],gust_params[i]['th_rad']*180/np.pi)
ax.set_ylabel('$\Theta_{rad} \ [^\circ]$')
ax.set_xlabel('$r/R$')
ax.legend(leg_lab)
ax.grid()
