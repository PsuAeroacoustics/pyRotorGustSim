import numpy as np
from scipy.fft import fft,ifft
from scipy.optimize import Bounds
import resonator as res

#%%

def get_sample_info(x,N_res,A_r,A_s,a_min,L_min):

    def get_res_geom(x):

        if N_res==1:
            B = np.array([1])
        
        else:
            # number of unique resonators per sample
            t = np.linspace(0,1, N_res)
            B = np.insert(x[2:],len(x[2:]),1)@np.array(((1-t)**2,2*t*(1-t),t**2))

        a = (x[0]-a_min)*B+a_min
        L = (x[1]-L_min)*B+L_min
        # a = x[0]*B
        # L = x[1]*B

        N = np.floor(A_r*A_s/np.sum(np.pi*a**2,axis = 0))

        return N,a,L

    def get_res_dist(x):

        N = len(A_s)*4
        t = np.arange(N)
        f_min = 1/(2*N)
        f_max = 1/8

        t2 = t/N
        B = x@np.array(((1-t2)**2,2*t2*(1-t2),t2**2))

        filt_ind = (0.5*np.sign(np.round(np.sin(2*np.pi*f_max*B*t),3))+0.5)[2::4]

        return filt_ind

    if len(x)==4:
        N,a,L = get_res_geom(x)

    elif len(x)==5:
        filt_ind = get_res_dist(x[2:])
        a = np.array([x[0]])
        L = np.array([x[1]])
        N = np.floor(A_r*A_s/(np.pi*a**2))
        N[filt_ind==0] = 0

    else:
        N,a,L = get_res_geom(x[:4])
        filt_ind = get_res_dist(x[4:])
        N[filt_ind==0] = 0

    return N,a,L

            
def init_res(f,a_n,L_n,a_c,L_c):
    '''
    This function creates and returns a resonator object with its complex impedance evaluated
    '''
    res_temp = res.resonator(a_n = a_n,L_n =L_n,a_c = a_c, L_c = L_c)
    res_temp.set_Z(f,model = 'k',rad = False,interior = False,loss = True,table = False)
    return res_temp

def init_fs(f,t_fs,r_fs,phi_fs,SPL,M,Z_cav):
    '''
    This function creates and returns a resonator object with its complex impedance evaluated
    '''
    fs_temp = res.fs(t = t_fs,r = r_fs,phi = phi_fs)
    fs_temp.set_Z(f,M = M,SPL = SPL,model = '2P',Z_cav = Z_cav)
    return fs_temp


def smeared_Z(f,N,a,L,N_res,facesheet,t_fs,phi_fs,N_fs,A_r,A_s,M,SPL):
    '''
    This function computes the smeared impedance and absorption for a sample that consists of multiple resonator cavities that have different geometries. 
    
    Parameters:
    f: frequency array [Hz]
    res_params: a nested array whose elements corresponding to each resonator in the sample. Each of the nested arrays must be of size (5,) and specifies
    [[# of this resonator in the sample, a_n,L_n,a_c,L_c]]

    Return:
    Z_tot: total non-dimensionalized (rho*c0) complex impedance of the sample (inverse of the total admittance)
    alpha: total absorption of the sample
    '''
    # res_params = np.array([[1.14238749e+02, 1.22225632e-03, 2.55143525e-01],[1.49012815e+02, 9.76120148e-04 ,4.66389696e-02],[3.12709586e+01, 1.38153591e-03, 1.25022033e-01],[5.91861319e+01, 1.05103046e-03, 1.85281711e-01]])

    # gets complex acoustic impedance of each resonator in a single sample/blade element
    # if N_res==1:
    #     Z = init_res(f,a_n = a,L_n = L/2,a_c = a,L_c = L/2).Z
    #     Z_tot = np.array([(N[i]*np.pi*a**2/A_s[i]*Z**-1)**-1 for i in range(elements)])

    # else:

    Z = np.array([init_res(f,a_n = a[i],L_n = L[i]/2,a_c = a[i],L_c = L[i]/2).Z for i in range(N_res)]).T
    Z_tot = (((N*np.pi*np.expand_dims(a,axis = -1)**2/A_s).T@Z.T**-1)**-1).T

    # if facesheet:

    #     # if linear:
    #     #     M = 0
    #     #     SPL = 60

    #     # else:
    #     #     SPL = 20*np.log10(np.sqrt(np.mean((z_loading_store[:,0]/A_s)**2,axis = 0))/20e-6)
    #     #     M = np.array(r)*R*omega/sos

    #     if N_res==1:
            
    #         r_fs = np.sqrt(phi_fs*np.pi*a**2/(N_fs*np.pi))
            
    #         if linear:
    #             Z_fs = init_fs(f,t_fs = t_fs,r_fs = r_fs,phi_fs = phi_fs,M = M,SPL = SPL,Z_cav = Z).Z
    #             Z_tot = (N*np.pi*a**2/A_s*np.expand_dims(Z+Z_fs,axis = -1)**-1).T**-1
    #         else:
    #             Z_fs = np.array([init_fs(f,t_fs = t_fs,r_fs = r_fs,phi_fs = phi_fs,M = M[i],SPL = SPL[i],Z_cav = Z).Z for i in range(elements)])
    #             Z_tot = (N*np.pi*a**2/A_s*(Z+Z_fs).T**-1).T**-1

    #     else:
    #         r_fs = np.min(np.sqrt(phi_fs*np.pi*a**2/(N_fs*np.pi)))
    #         N_fs = np.floor(phi_fs*np.pi*a**2/(np.pi*r_fs**2))
    #         phi_fs = N_fs*np.pi*r_fs**2/(np.pi*a**2)

    #         # r_fs = np.sqrt((np.min(np.pi*a**2)*phi)/np.pi)
    #         # N_fs = np.floor(phi*a**2/r_fs**2)
    #         # phi = N_fs*np.pi*r_fs**2/(np.pi*a**2)
    #         if linear:
    #             Z_fs = np.array([init_fs(f,t_fs = t_fs,r_fs = r_fs,phi_fs = phi_fs[i],M = M,SPL = SPL,Z_cav = Z[i]).Z for i in range(N_res)]).T
    #         else:
    #             Z_fs = np.array([[init_fs(f,t_fs = t_fs,r_fs = r_fs,phi_fs = phi_fs[i],M = M[ii],SPL = SPL[ii],Z_cav = Z[i]).Z for i in range(N_res)] for ii in range(elements)]).transpose(0,-1,1)


    alpha = 1 - abs((Z_tot-1)/(Z_tot+1))**2


    return Z_tot, alpha

def get_filt_resp(Z_tot,odd = True):

    filt_resp_shape = Z_tot.shape
    # reflection coefficient
    R = (Z_tot-1)/(Z_tot+1)
    R[Z_tot==1] = 1+1j*0
    
    if odd:
        filt_resp = np.ones((filt_resp_shape[0]*2+1,filt_resp_shape[-1]),dtype = complex)
        filt_resp[1:int(filt_resp_shape[0])+1] = R
        filt_resp[int(filt_resp_shape[0])+1:] = np.conj(R)[::-1]

    else:
        # complex admittance of the resonators
        filt_resp = np.ones((filt_resp_shape[0]*2,filt_resp_shape[-1]),dtype = complex)
        filt_resp[1:int(filt_resp_shape[0])] = R[:-1]
        filt_resp[int(filt_resp_shape[0])+1:] = np.conj(R[:-1])[::-1]

    return filt_resp

def apply_filt(loads,filt_resp):

    Xm = fft(loads,axis = 0)
    # Xm_filt = np.zeros(Xm.shape,dtype=complex)
    # Xm_filt = (Xm.T*filt_resp).T
    Xm_filt = (Xm*np.expand_dims(filt_resp,axis = -1))
    # performs an ifft to convert filtered loads back to the time domain
    xn_filt = np.real(ifft(Xm_filt,axis =0))
    return xn_filt

def get_opt_bounds(x,r_min,r_max,L_min,L_max):
    if len(x) == 4:
        bounds = Bounds(lb = (r_min,L_min,0,0), ub = (r_max,L_max,1,1),keep_feasible=True)
    elif len(x) == 5:
        bounds = Bounds(lb = (r_min,L_min,-1,-1,-1), ub = (r_max,L_max,1,1,1),keep_feasible=True)
    else:
        bounds = Bounds(lb = (r_min,L_min,0,0,-1,-1,-1), ub = (r_max,L_max,1,1,1,1,1),keep_feasible=True)
    return bounds
 
