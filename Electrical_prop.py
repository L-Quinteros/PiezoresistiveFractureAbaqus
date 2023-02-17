# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 21:32:50 2021
@author: Enrique GM
"""
import numpy as np
import math

def Eff_conductividy(dco,Lambdao,L_CNT,d_CNTo,sigma_iter,sigma_M,vio):

    # Properties of CNTs
    # ==================
    # Length
    L_CNT = L_CNT*10**(-6);
    # Diameter
    d_CNT = d_CNTo*10**(-9);
    # Conductivity of CNTs (L=longitudinal, T=transversal)
    sigmaL_CNT = sigma_iter;
    sigmaT_CNT = sigma_iter;
    CNT_prop = [d_CNT,L_CNT];
    
    # MICROMECHANICS MODELLING
    # ==========================
    # PERCOLATION THRESHOLD
    s = L_CNT/d_CNT;
    I = 1.0/(1.27327);
    fc = np.pi/(5.77*s*I);
    # Volume fraction
    vi = vio;  # Volume fraction
    # 1. INTERPHASE
    [Sigma_int_EH,t_EH,Sigma_int_CN,t_CN] = interphase_CNT(CNT_prop,vi,fc,dco,Lambdao);
    Sigma_int = np.array([Sigma_int_EH,Sigma_int_CN]);
    t = np.array([t_EH,t_CN]);
    # 2. EQUIVALENT CYLINDER
    [sigmaT_EH,sigmaL_EH,sigmaT_CN,sigmaL_CN] = equivalent_filler(sigmaT_CNT,sigmaL_CNT,Sigma_int,L_CNT,t,d_CNT);
    equi_filler = np.array([sigmaT_EH,sigmaL_EH,sigmaT_CN,sigmaL_CN])
    # 3. MICROMECHANICS PIEZOELECTRIC CNT
    [sigma_EFF,Xi] = micro_piezoCNT_cf(equi_filler,vi,CNT_prop,sigma_M,fc,t,np.array([0,0,0]));
    return sigma_EFF


def interphase_CNT(CNT_prop,f,fc,dco,Lambdao):

    d_CNT=CNT_prop[0];
    L_CNT=CNT_prop[1];
    rc=d_CNT/2.;
    
    # Assumptions
    # =================
    # Maximum possible distance that allows the tunneling penetration of
    # electrons
    dc = dco*10**(-9);
    # Planck constant
    hplanck = 6.626068*10**(-34);
    # Mass of an electron
    m = 9.10938291*10**(-31);
    #% Electric charge of an electron
    ee = 1.602176565*10**(-19);
    # Potential barrier height
    Lambda=Lambdao*ee;
    
    # Initial calculations
    # ==================
    # Aspect ratio
    asp = L_CNT/d_CNT;
    # Contact area of CNTs
    a_EH=np.pi*(d_CNT/2)**2;
    a_CN=np.pi*(d_CNT/2)**2;    
    
    # Percentage of percolated CNTs *)
    Xi = (f**(1./3.)-fc**(1./3.))/(1-fc**(1./3.));
    
    
    # Interphase
    # ==================
    # CONDUCTIVE NETWORK:
    da=(fc/f)**(1./3.)*dc;
    # Tunneling-type contact resistance between two CNTs
    Rint_CN = (da*hplanck**2./(a_CN*ee**2.*(2.*m*Lambda)**0.5))*np.exp((4.*np.pi*da/hplanck)*(2.*m*Lambda)**0.5);
    # Thickness of the interphase *)
    t_CN = da/2.;    # Power law relation after percolation
    # Conductivity of the interphase
    Sigma_int_CN = da/(a_CN*Rint_CN);
    
    # ELECTRON HOPPING:
    # Tunneling-type contact resistance between two CNTs
    Rint_EH = (dc*hplanck**2./(a_EH*ee**2.*(2.*m*Lambda)**0.5))*np.exp((4.*np.pi*dc/hplanck)*(2.*m*Lambda)**0.5);
    # Thickness of the interphase *)
    t_EH = dc/2.;    # Constant distance before percolation
    # Conductivity of the interphase
    Sigma_int_EH = dc/(a_EH*Rint_EH);
    
    
    if Xi>=0:
        Sigma_int_EH = Sigma_int_CN;
        t_EH = t_CN;

    
    
    return [Sigma_int_EH,t_EH,Sigma_int_CN,t_CN]



def equivalent_filler(sigmaT_CNT,sigmaL_CNT,Sigma_int,L,t,d_CNT):

    # Radio of the CNT
    rc=d_CNT/2;
    
    # Electron hopping
    # ================
    sigmaT_EH=(Sigma_int[0]/(L+2.*t[0]))*(L*(2.*rc**2*sigmaT_CNT+(sigmaT_CNT+Sigma_int[0])*(t[0]**2.+2.*rc*t[0]))/(2*rc**2*Sigma_int[0]+(sigmaT_CNT+Sigma_int[0])*(t[0]**2.+2.*rc*t[0]))+2.*t[0]);
    sigmaL_EH=(L+2.*t[0])*Sigma_int[0]*(sigmaL_CNT*rc**2.+Sigma_int[0]*(2.*rc*t[0]+t[0]**2.))/(2.*sigmaL_CNT*rc**2.*t[0]+2.*Sigma_int[0]*(2.*rc*t[0]+t[0]**2)*t[0]+Sigma_int[0]*L*(rc+t[0])**2.);
         
    # Conductive Networks
    # ================
    sigmaT_CN=(Sigma_int[1]/(L+2.*t[1]))*(L*(2.*rc**2*sigmaT_CNT+(sigmaT_CNT+Sigma_int[1])*(t[1]**2.+2.*rc*t[1]))/(2*rc**2*Sigma_int[1]+(sigmaT_CNT+Sigma_int[1])*(t[1]**2.+2.*rc*t[1]))+2.*t[1]);
    sigmaL_CN=(L+2.*t[1])*Sigma_int[1]*(sigmaL_CNT*rc**2.+Sigma_int[1]*(2.*rc*t[1]+t[1]**2.))/(2.*sigmaL_CNT*rc**2.*t[1]+2.*Sigma_int[1]*(2.*rc*t[1]+t[1]**2)*t[1]+Sigma_int[1]*L*(rc+t[1])**2.);
     
    return [sigmaT_EH,sigmaL_EH,sigmaT_CN,sigmaL_CN]



def micro_piezoCNT_cf(equi_filler,vi,CNT_prop,sigma_M,fc,t,strain):


    # CNT
    # ================================
    d=CNT_prop[0];
    L=CNT_prop[1];
    rc=d/2.;
    
    # Approach of the equivalent fiber
    # ================================
    sigmaT_EH = equi_filler[0];
    sigmaL_EH = equi_filler[1];
    sigmaT_CN = equi_filler[2];
    sigmaL_CN = equi_filler[3];
    
    # Percentage of percolated CNTs
    if vi<fc:
        Xi=0.;
    else:
        Xi = (vi**(1./3.)-fc**(1./3.))/(1.-fc**(1./3.));
    
    # ================================
    # EFFECTIVE PROPERTIES
    #  ================================
    
    # Effective volume fraction
    feff_EH = vi*(rc+t[0])**2*(L+2*t[0])/(rc**2*L);
    if Xi>=0:
        feff_CN = vi*(rc+t[1])**2.*(L+2.*t[1])/(rc**2.*L);
    
    # Construction of the conductivity tensors
    # Electrical conductivity tensor of the effective filler: Electron hopping *)
    Sigma_EH = np.array([[sigmaL_EH, 0., 0.],
        [0.,sigmaT_EH, 0.],
        [0.,0.,sigmaT_EH]]);
    # Electrical conductivity tensor of the effective filler: Conductive network
    Sigma_CN = np.array([[sigmaL_CN, 0., 0.],
        [0.,sigmaT_CN, 0.],
        [0.,0.,sigmaT_CN]]);
    # Electrical conductivity tensor of the matrix
    Sigma_M = np.array([[sigma_M, 0., 0.],
        [0.,sigma_M, 0.],
        [0.,0.,sigma_M]]);
    
    # Eshelby's tensor
    # ================================
    # Aspect ratio of the equivalent filler: Electron hopping
    Are = (L + 2.*t[0])/(2.*rc + 2.*t[0]);
    S22 = (Are)*(Are*(Are**2. - 1.)**0.5 - math.acosh(Are))/(2.*(Are**2. - 1.)**(3./2.));
    S33 = S22;
    S11 = 1. - 2.*S22;
    SEH = np.array([[S11, 0., 0.],
        [0., S22, 0.],
        [0., 0., S33]]);
    # Aspect ratio of the equivalent filler: Conductive Network
    SCN =np.array([[0., 0., 0.],
        [0.,0.5,0.],
        [0.,0,0.5]]);
    
    # Pre-Rotation of matrices
    # ================================
    Sigma_EH=((rot(0.,np.pi/2.,0.)*Sigma_EH)*rot(0.,np.pi/2.,0.).T)
    SEH=(rot(0.,np.pi/2.,0.)*SEH)*rot(0.,np.pi/2.,0.).T
    Sigma_CN=((rot(0.,np.pi/2.,0.)*Sigma_CN)*rot(0.,np.pi/2.,0.).T);
    SCN=((rot(0.,np.pi/2.,0.)*SCN)*rot(0.,np.pi/2.,0.).T);
    
    
    
    # ELECTRON HOPPING
    # ================================
    # Field concentration factor
    delta = np.eye(3,3);
    TEH = (delta+(SEH*np.linalg.inv(Sigma_M))*(Sigma_EH-Sigma_M));
    AdilEH = np.linalg.inv(TEH);
    AdilEHoa = np.real(Orientational_average_closed_form(AdilEH,strain));
    
    # CONDUCTIVE NETWORKS
    # ================================
    # Field concentration factor
    T = (delta+(SCN*np.linalg.inv(Sigma_M))*(Sigma_CN-Sigma_M));
    TCN = T;
    AdilCN = np.linalg.inv(TCN);
    AdilCNa = np.real(Orientational_average_closed_form(AdilCN,strain));
    
    
    # ================================
    # EFFECTIVE PROPERTIES
    # ================================
    EHM = feff_EH*Orientational_average_closed_form((Sigma_EH-Sigma_M)*AdilEH*np.linalg.inv(((1.-feff_EH)*delta+feff_EH*AdilEHoa)),strain);
    CNM = feff_CN*Orientational_average_closed_form((Sigma_CN-Sigma_M)*AdilCN*np.linalg.inv(((1.-feff_CN)*delta+feff_CN*AdilCNa)),strain);
    Sigma_eff = Sigma_M + (1.-Xi)*EHM+Xi*CNM;
    #print(AdilCNa)
    #print(((1.-feff_CN)*delta+feff_CN*AdilCNa))
    #print(((1.-feff_CN)*delta+feff_CN*AdilCNa))
    #print(np.linalg.inv(((1.-feff_CN)*delta+feff_CN*AdilCNa)))
    
    
    return [np.real(Sigma_eff),Xi]


def rot(Beta,Alpha,Psi):

    # Compute Rotation tensor
    #
    # INPUT:
    # Beta -> Rotation around x1 axis
    # Alpha -> Rotation around x2 axis
    #
    # OUTPUT:
    # Q -> Rotation matrix of 3x3 tensor
    # Author: E. García-Macias
    # ------------------------------------------------------------------------
    
    orden=np.array([0,2,1]);
    
    R = np.zeros((3,3,3))
    
    R[0,:,:] = np.array([[1.,0.,0.],[0.,np.cos(Beta),np.sin(Beta)],[0.,-np.sin(Beta),np.cos(Beta)]]);
    R[1,:,:] = np.array([[np.cos(Alpha),0.,-np.sin(Alpha)],[0.,1.,0.],[np.sin(Alpha),0.,np.cos(Alpha)]]);
    R[2,:,:] = np.array([[np.cos(Psi),np.sin(Psi),0.],[-np.sin(Psi),np.cos(Psi),0.],[0.,0.,1.]]);
    
    Q = (np.mat(np.squeeze(R[orden[2],:,:]))*np.mat(np.squeeze(R[orden[1],:,:])))*np.mat(np.squeeze(R[orden[0],:,:]))
    
    return  Q



def Orientational_average_closed_form(Sigma,strain):

    Sigma11 = Sigma[0,0];
    Sigma33 = Sigma[2,2];
    Eps1 = strain[0]+1.
    Eps2 = strain[1]+1.
    Eps3 = strain[2]+1.
    
    SigmaPromedio11=(1./15.)*(5.*(2.*Sigma11+Sigma33)+4.*5.**(1./2.)*np.pi**2.*(\
      Sigma11+(-1.)*Sigma33)*((1j*(1./8008.))*5.**(-1./2.)*\
      Eps1*(Eps1+(-1.)*Eps2)*Eps2*(Eps1**2.*Eps2**2.)**(1./2.)*\
     Eps3**2.*(5.*(65760.+(-37329.)*Eps2+Eps1*((-37329.)+17192.*\
      Eps2))+(-40.)*(11448.+(-6521.)*Eps2+Eps1*((-6521.)+3024.*\
      Eps2))*Eps3+8.*(21860.+(-12474.)*Eps2+27.*Eps1*((-462.)+\
      215.*Eps2))*Eps3**2)*np.pi**(-1.)*(np.log((1j*(-1.))*\
      Eps1*Eps2**(-1.))+(-1.)*np.log(1j*Eps1*Eps2**(-1.)))**(\
      -1.)+(1j*(-1./8008.))*5.**(-1./2.)*Eps1*Eps2*(Eps1**2.*\
      Eps2**2.)**(1./2.)*Eps3**2.*(5.*Eps2*((-21920.)+12443.*Eps2)+\
      40.*(5480.+(-3816.)*Eps2+741.*Eps2**2.)*Eps3+(-2.)*(62215.+8.*\
      Eps2*((-7225.)+2079.*Eps2))*Eps3**2.+40.*Eps1*((-2740)+\
      7632.*Eps2+(-3631.)*Eps2**2+72.*((-53.)+14.*Eps2**2.)*Eps3+(\
      2890+9.*Eps2*((-224)+43.*Eps2))*Eps3**2)+Eps1**2.*(62215.+\
      29640.*Eps3+8.*(Eps2*((-18155.)+8316.*Eps2)+90.*(56.+(-43.)*\
      Eps2)*Eps2*Eps3+9.*((-462.)+215.*Eps2)*Eps3**2.)))*np.pi**(\
      -1.)*(np.log((1j*(-1.))*Eps1*Eps2**(-1.))+(-1.)*np.log(1j\
      *Eps1*Eps2**(-1.)))**(-1.)));
    
    
    SigmaPromedio22=(1./15.)*(5.*(2.*Sigma11+Sigma33)+4.*5.**(1./2.)*np.pi**2.*(\
      Sigma11+(-1)*Sigma33)*((1j*(-1./8008.))*5.**(-1./2.)*\
      Eps1*(Eps1+(-1.)*Eps2)*Eps2*(Eps1**2.*Eps2**2.)**(1./2.)*\
     Eps3**2.*(5.*(65760+(-37329)*Eps2+Eps1*((-37329)+17192.*\
      Eps2))+(-40.)*(11448+(-6521)*Eps2+Eps1*((-6521)+3024.*\
      Eps2))*Eps3+8.*(21860+(-12474)*Eps2+27.*Eps1*((-462)+\
      215.*Eps2))*Eps3**2.)*np.pi**(-1.)*(np.log((1j*(-1))*\
      Eps1*Eps2**(-1.))+(-1)*np.log(1j*Eps1*Eps2**(-1)))**(\
      -1)+(1j*(-1./8008.))*5.**(-1./2.)*Eps1*Eps2*(Eps1**2.*\
      Eps2**2)**(1./2.)*Eps3**2.*(5.*Eps2*((-21920)+12443.*Eps2)+\
      40.*(5480+(-3816)*Eps2+741.*Eps2**2.)*Eps3+(-2)*(62215+8.*\
      Eps2*((-7225)+2079.*Eps2))*Eps3**2+40.*Eps1*((-2740)+\
      7632.*Eps2+(-3631)*Eps2**2.+72.*((-53)+14.*Eps2**2)*Eps3+(\
      2890+9.*Eps2*((-224)+43.*Eps2))*Eps3**2)+Eps1**2.*(62215+\
      29640.*Eps3+8.*(Eps2*((-18155)+8316.*Eps2)+90.*(56+(-43)*\
      Eps2)*Eps2*Eps3+9.*((-462)+215.*Eps2)*Eps3**2)))*np.pi**(\
      -1)*(np.log((1j*(-1))*Eps1*Eps2**(-1))+(-1)*np.log(1j\
      *Eps1*Eps2**(-1)))**(-1)));
    
    
    SigmaPromedio33=(1./15.)*(5.*(2.*Sigma11+Sigma33)+(1j*(-1./1001.))*\
      Eps1*Eps2*(Eps1**2.*Eps2**2)**(1./2.)*Eps3**2.*(5.*Eps2*((\
      -21920)+12443.*Eps2)+40.*(5480+(-3816)*Eps2+741.*Eps2**2)*\
      Eps3+(-2)*(62215+8.*Eps2*((-7225)+2079.*Eps2))*Eps3**2+\
      40.*Eps1*((-2740)+7632.*Eps2+(-3631)*Eps2**2+72.*((-53)+\
      14.*Eps2**2)*Eps3+(2890+9.*Eps2*((-224)+43.*Eps2))*\
     Eps3**2)+Eps1**2.*(62215+29640.*Eps3+8.*(Eps2*((-18155)+\
      8316.*Eps2)+90.*(56+(-43)*Eps2)*Eps2*Eps3+9.*((-462)+\
      215.*Eps2)*Eps3**2)))*np.pi*((-1)*Sigma11+Sigma33)*(np.log((\
      1j*(-1))*Eps1*Eps2**(-1))+(-1)*np.log(1j*Eps1*\
      Eps2**(-1)))**(-1));
    
    sigmaeff = np.array([[SigmaPromedio11,0.,0.],[0.,SigmaPromedio22,0.],[0.,0.,SigmaPromedio33]])
    
    return sigmaeff 


def Piezoresistivity(dco,Lambdao,L_CNT,d_CNT,sigma_iter,sigma_M,vio,str_comp,str_tens):

    
    # Properties of CNTs
    # ==================
    # Length
    L_CNT = L_CNT*10**(-6);
    # Diameter
    d_CNT = d_CNT*10**(-9);
    # Conductivity of CNTs (L=longitudinal, T=transversal)
    sigmaL_CNT = sigma_iter;
    sigmaT_CNT = sigma_iter;
    CNT_prop = [d_CNT,L_CNT];
    # Volume fraction
    vio = vio;  # Transformation to volume fraction
    
    # Range of deformations
    strain_vector = np.array([np.linspace(float(str_comp),float(str_tens),699, endpoint=True)/100.]).T;
    strain_vector = np.row_stack((np.array(0),strain_vector))
    
    # Initialize variables
    Drho_12 = np.zeros((len(strain_vector),1));
    Drho_11 = Drho_12;
    Xi = np.zeros((1,len(strain_vector)));
    fc = np.zeros((1,len(strain_vector)));
    
    for j in np.arange(0,len(strain_vector)):
        
        # Uni-axial stretching
        # ====================
        strain = np.array([0.,0.,float(strain_vector[j])]); # Constrained lateral displacement
        strainvol = 1.+strain;
        # Stretching induced volume expansion
        vi = float(vio/(strainvol[0]*strainvol[1]*strainvol[2]))
        
        # PERCOLATION THRESHOLD
        s = L_CNT/d_CNT;
        I = 1./(1.27327+0.25457*strainvol[2]-0.25461*np.log(strainvol[2]+1.));
        fc[0,j] = np.pi/(5.77*s*I);
        # 1. INTERPHASE
        [Sigma_int_EH,t_EH,Sigma_int_CN,t_CN] = interphase_CNT(CNT_prop,vi,fc[0,j],dco,Lambdao);
        Sigma_int = [Sigma_int_EH,Sigma_int_CN];
        t = np.array([t_EH,t_CN])
        # 2. EQUIVALENT CYLINDER
        [sigmaT_EH,sigmaL_EH,sigmaT_CN,sigmaL_CN] = equivalent_filler(sigmaT_CNT,sigmaL_CNT,Sigma_int,L_CNT,t,d_CNT);
        equi_filler = np.array([sigmaT_EH,sigmaL_EH,sigmaT_CN,sigmaL_CN])
        # 3. MICROMECHANICS PIEZORESISTIVE CNT
        [sigma_EFF,Xi[0,j]] = micro_piezoCNT_cf(equi_filler,vi,CNT_prop,sigma_M,fc[0,j],t,strain);
        
        # Xi/Xi0
        if j==0:
            sigma_EFFo=sigma_EFF[1,1];
        Drho_12[j] = (sigma_EFFo/sigma_EFF[0,0])-1;
        Drho_11[j] = (sigma_EFFo/sigma_EFF[2,2])-1;
        
    
    # Characterization of the strain-sensing curves
    
    # L11
    # ---------------------------------------
    sensitivity = Drho_11[1:].T;
    strainserie = strain_vector[1:].T
    # Compression
    poscompress = np.argwhere(strainserie>=0)-1;
    x = -(strainserie[0,1:poscompress[0,1]])
    x=x[::-1]
    y=-(sensitivity[0,1:poscompress[0,1]])
    
    y=y[::-1]
    #x = -np.flip(strainserie[0,1:poscompress[0,1]]);
    #y = -np.flip(sensitivity[0,1:poscompress[0,1]]);
    
    L11_comp = gauge_reg(x,y);
    # Traction

    x = strainserie[0,poscompress[0,1]+1:];
    y = sensitivity[0,poscompress[0,1]+1:];
    
    L11_tract = gauge_reg(x,y);
            
    # L12
    # ---------------------------------------
    sensitivity = Drho_12[1:].T
    strainserie = strain_vector[1:].T
    # Compression
    poscompress = np.argwhere(strainserie>=0)-1;
    
    x = -(strainserie[0,0:poscompress[0,1]])
    x=x[::-1]
    y=-(sensitivity[0,0:poscompress[0,1]])
    y=y[::-1]
    
    #x = -np.flip(strainserie[0,0:poscompress[0,1]]);
    #y = -np.flip(sensitivity[0,0:poscompress[0,1]]);
    L12_comp = gauge_reg(x,y);
    # Traction
    x = strainserie[0,poscompress[0,1]+1:-1];
    y = sensitivity[0,poscompress[0,1]+1:-1];
    L12_tract = gauge_reg(x,y);
    
    # L44
    L44_tract = (L11_tract-L12_tract)/2;
    L44_comp = (L11_comp-L12_comp)/2; 
    
    
    # Output
    Xi = np.delete(Xi,0)
    fc = np.delete(fc,0)
    Drho_12 = np.delete(Drho_12,0)
    Drho_11 = np.delete(Drho_11,0)
    strain_vector = np.delete(strain_vector,0)

    return [strain_vector,Drho_11,Drho_12,L11_tract,L12_tract,L11_comp,L12_comp,L44_tract,L44_comp,fc,Xi]




def gauge_reg(x,y):
    
    if x[0]!=0:
      #xtot = np.pad(x, (1,0), 'constant')
      #ytot = np.pad(y, (1,0), 'constant')
      xtot=np.zeros([len(x)+1])
      ytot=np.zeros([len(x)+1])
      xtot=np.zeros([len(x)+1])
      ytot=np.zeros([len(y)+1])
      xtot[0]=0
      ytot[0]=0
      xtot[1:]=x
      ytot[1:]=y
    else:
      xtot = x;
      ytot = y; 
    
    p = np.polyfit(xtot[0:2],ytot[0:2],1);
    yfit = np.polyval(p,xtot[0:2]);
    yresid = ytot[0:2]-yfit;
    SSresid = sum(yresid**2.);
    
    SStotal = (len(ytot)-1.)*np.var(ytot);
    
    rsq = 1.-SSresid/SStotal;
    
    rsq_limit = 0.9999;
    cont = 1;
    while rsq>rsq_limit and float(cont+1)<=len(xtot):
      cont=cont+1;
        
      x = xtot[0:cont]
      y = ytot[0:cont]
    
      p = np.polyfit(x,y,1);
      yfit = np.polyval(p,x);
      yresid = y-yfit;
      SSresid = sum(yresid**2);
      SStotal = (len(y)-1)*np.var(y);
      rsq = 1-SSresid/SStotal;

    
    p_ref = ytot[cont-1]/(xtot[cont-1]);
    
    
    
    cont_level = 0;
    cont_level2 = 0;
    length1 = xtot[-1];
    length5 = xtot[-1];
    for i in np.arange(cont,len(xtot)):
    
        if np.abs(1-ytot[i]/(xtot[i]*p_ref))>0.05 and cont_level==0:
          length5 = xtot[i]
          cont_level = 1;
        if abs(1-ytot[i]/(xtot[i]*p_ref))>0.01 and cont_level2==0:
          length1 = xtot[i]
          cont_level2 = 1;
    
    
    # Store variables
    
    gauge = p_ref;
    linlength1 = length1;
    linlength5 = length5;
    return gauge








def electrical_prop(L_cnt,D_cnt,dco,Lambdao,sigma_m,sigma_cnt,vf,str_comp,str_tens):
    L_cnt=L_cnt*1e6
    D_cnt=D_cnt*1e9
    dco=dco*1e9
    sigmaEFF = Eff_conductividy(dco,Lambdao,L_cnt,D_cnt,sigma_cnt,sigma_m,vf)
    [strain_vector,Drho_11,Drho_12,L11_tract,L12_tract,L11_comp,L12_comp,L44_tract,L44_comp,fc,Xi] = Piezoresistivity(dco,Lambdao,L_cnt,D_cnt,sigma_cnt,sigma_m,vf,str_comp,str_tens)
    return sigmaEFF,L11_tract,L12_tract,L44_tract










# INPUT PARAMETERS
# Matrix phase [Isotropic material]
#densm = 1.12               # Density g/cm3
sigma_m = 1.036000000000000e-10        # Electrical conductivity S/m
 
# MWCNTs [Transverse isotropic material]
# Mechanical properties
#densp = 1.392586934259224                # Density g/cm3 
L_cnt = 3.21e-6               # Length [microns] 
D_cnt = 10.35e-9              # Diameter [nm] 
# Electrical properties
dco = 0.22e-9     # Interparticle distance [nm]
Lambdao = 0.69 # Height of the potential barrier [eV]
# Electrical conductivities of CNTs to be analysed
sigma_cnt = 10**2
    
# ELECTRICAL PROPERTIES
#vf=0.01
# STRAIN SENSING CURVES
# Maximum compression
str_comp = -5;  # [%]
# Maximum traction
str_tens = 5;  # [%]

#sigmaEFF,L11_tract,L12_tract,L44_tract=electrical_prop(L_cnt,D_cnt,dco,Lambdao,sigma_m,sigma_cnt,vf,str_comp,str_tens)



AspectRatio = np.array([100,200,300,400,500])
VolFrac = np.linspace(0.00000001,0.05,1001)

ValEff = np.zeros([1001,5])
L11 = np.zeros([1001,5])


for i, AR in enumerate(AspectRatio):
    for j, vf in enumerate(VolFrac):
        print(i)
        L_cnt = AR*D_cnt
        sigmaEFF,L11_tract,L12_tract,L44_tract=electrical_prop(L_cnt,
                                                               D_cnt,dco,
                                                               Lambdao,sigma_m,
                                                               sigma_cnt,vf,
                                                               str_comp,
                                                               str_tens)
        L11[j,i] = L11_tract
        ValEff[j,i] = sigmaEFF[0,0]
    
    
    


