#Fracture energy calculation

import numpy as np 



def L_ctheta(D_cnt,sigma_ult,tau_int,theta,Ac,mu):
    num=D_cnt*(sigma_ult*(1-Ac*np.tan(theta)))
    den=tau_int*np.exp(mu*theta)*4
    return num/den    

def g_dis(theta,p,q,theta_min,theta_max):
    Nint=101
    num=lambda theta: (np.sin(theta)**(2*p-1.))*(np.cos(theta)**(2*q-1.))
    x=np.linspace(theta_min,theta_max,Nint) 
    int1=np.zeros([Nint])
    for i in range(len(x)):
        int1[i]=num(x[i])
    den = np.trapz(int1,x)    
        
    return num(theta)/den


def Fracture_potential(l,D_cnt,sigma_ult,Ac,theta,tau_int,mu,L_cnt,E_cnt):
    L_ct=L_ctheta(D_cnt,sigma_ult,tau_int,theta,Ac,mu)
    w=0
    if l<L_ct:
        w=(l**2)*tau_int*np.pi*D_cnt*np.exp(mu*theta)/2
    else:
        w=np.pi*(D_cnt**2)*(sigma_ult**2)*L_cnt/(8*E_cnt)
    return w

def integrand(theta,l):
    func = lambda theta: np.cos(theta)*g_dis(theta,p,q,theta_min,theta_max)*((8*f_p)/(np.pi*(D_cnt**2)*L_cnt))
    w= Fracture_potential(l,D_cnt,sigma_ult,Ac,theta,tau_int,mu,L_cnt,E_cnt)*func(theta)
    #w= Fracture_potential(l,D_cnt,sigma_ult,Ac,theta,tau_int,mu,L_cnt,E_cnt)
    return w
    
def Energy_release_rate(theta_min,theta_max,L_min,L_max,f_p,D_cnt,L_cnt,sigma_ult,Ac,tau_int,mu,E_cnt,G_m,p,q):
    Nint=200
    theta_serie=np.linspace(theta_min,theta_max,Nint)
    L_serie = np.linspace(L_min,L_max,Nint)
    int2=np.zeros([len(theta_serie)])
    for i in range(len(theta_serie)): 
       theta=theta_serie[i]   
       #Initialize variable
       int1=np.zeros([len(L_serie)])
       for j in range(len(L_serie)):
           l=L_serie[j]
           factor=((8*f_p*np.cos(theta))/((D_cnt**2)*np.pi*L_cnt))
           g_f=g_dis(theta,p,q,theta_min,theta_max)   
           int1[j]=factor*Fracture_potential(l,D_cnt,sigma_ult,Ac,theta,tau_int,mu,L_cnt,E_cnt)*g_f   
       int2[i]=np.trapz(int1,L_serie)    
       
    return np.trapz(int2,theta_serie)+G_m
    

    
