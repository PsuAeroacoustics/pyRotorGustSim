import numpy as np

#%%

class Atmosphere:
    def __init__(self,rho = 1.125,sos = 340,nu =14.88e-6):
        self.rho = rho
        self.sos = sos
        self.nu = nu
        self.mu = self.nu/self.rho

class aircraft:
    def __init__(self,N_rotor = 1):
        self.N_rotor = N_rotor

class rotor:

    def __init__(self,Nb,R,e,c,th0,th_tw,N_elements,af,Cl_a,origin,omega,V_c,C_T_target,atmos):

        self.Nb = Nb
        self.R = R
        self.c = c
        self.e = e
        self.th0 = th0
        self.th_tw = th_tw
        self.origin = origin
        self.N_elements = N_elements
        self.sigma = self.Nb*self.c/(np.pi*self.R)

        elems = (np.arange(self.N_elements+1)*(self.R-self.e)/(self.N_elements)+self.e)/self.R
        r = 0.5*(elems[1:]+elems[:-1])
        c = self.c/self.R
        th = self.th0+r*self.th_tw

        self.omega = omega
        self.V_c = V_c
        self.C_T_target= C_T_target
        self.lam_c = V_c/(self.omega*self.R)

        Re = omega*R*r*c/atmos.nu
        M = omega*R*r/atmos.sos

        self.blades = [blade(r,self.c,self.R,th,af,Cl_a=Cl_a,offset = b_iter*2*np.pi/self.Nb,Re = Re,M = M) for b_iter in range(self.Nb)]

    def set_loads(self):
        self.CT =np.sum(np.array([np.trapz(b.dCT,x = b.r) for b in self.blades])) 
        self.CP =np.sum(np.array([np.trapz(b.dCP,x = b.r) for b in self.blades])) 


class blade:
    def __init__(self,r,c,R,th,af,Cl_a,offset,Re,M):
        self.r = r
        self.c = c
        self.R = R
        self.th = th
        self.af = af
        self.offset = offset
        self.Re = Re
        self.M = M
        self.Cl_a = Cl_a

    def set_twist(self,rotor):
        self.th = rotor.th0+self.r*rotor.th_tw

    def set_loads(self):

        self.U = np.sqrt(self.lam**2+self.r**2)
        self.phi = np.arctan2(self.lam,self.r)
        self.aoa = self.th-self.phi

        self.CL,self.CD = get_af_coeffs(self.af, self.aoa*180/np.pi, self.Re, self.M)

        dCz = self.CL*np.cos(self.phi)-self.CD*np.sin(self.phi)
        dCx = self.CL*np.sin(self.phi)+self.CD*np.cos(self.phi)

        self.dCT = 0.5*self.c/(np.pi*self.R)*self.U**2*dCz
        self.dCP = 0.5*self.c/(np.pi*self.R)*self.r*self.U**2*dCx

def get_af_coeffs(af, aoa, Re, M):
    aero = af.get_aero_from_neuralfoil(alpha=aoa.flatten(), Re=Re.flatten(),mach =M.flatten(),model_size='xlarge')
    CL,CD = aero['CL'].reshape(aoa.shape),aero['CD'].reshape(aoa.shape)
    return CL,CD

    #     return err
