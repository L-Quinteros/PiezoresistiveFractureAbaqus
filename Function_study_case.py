# -*- coding: utf-8 -*-
# Model electrically conductive strip
#######################################
# Import modules and libraries
#######################################
# Import Abaqus modules
from abaqus import *
from abaqusConstants import *
from caeModules import *
from abaqus import getInput
from odbAccess import *

##########################
import numpy as np
import sys
import fileinput
import math

##########################
# Import Study cases
from Studies import *

##########################
# Import post process
from Post_processing import post_process

########################################
# Import tool-box for electromechanical properties
from Electrical_prop import *
from Fracture_energy import *
from Mechanical_prop import *
import os

# Change the directory
os.chdir(
    r"/rds/general/user/ldq20/home/case_study_length/L_scale_study"
)

# Run the case



def run_the_case(vf,theta,Inc,case_number,Vimp,Dimp,file_name,k_value,n_value,meshsize,meshsize2):


    ###########################
    # RUN STUDY CASES
    ###########################
    ###########################
    ###ELASTIC TENSOR###
    ###########################
    
    # Geometrical properties
    ######################
    L_cnt = 3.20995854347668e-6  # Length [microns]
    D_cnt = 10.3538294770233e-9  # Diameter [nm]
    mu = 0  #
    Ac = 0.083  # Orientation limit angle
    # Matrix properties
    ######################
    E_matrix = 2.5e9  # Youngs modulus MPa
    nu_matrix = 0.28  # Poissons ratio
    sigma_m = 1.036000000000000e-10  # Electrical conductivity S/m
    G_m = 133  # [J/m^2]
    # CNT properties
    ######################
    E_cnt = 700e9  # Youngs modulus MPa
    nu_cnt = 0.3  # Poissons ratio
    sigma_ult = 35e9  # CNT ultimate strength [Pa]
    tau_int = 47e6  # CNT interfatial strength [Pa] (epoxy-CNT)
    sigma_cnt = 10 ** 2
    # Interphase - just for elastic modulus and poisson ratio
    ######################
    E_i = 2.17e9
    nu_i = nu_matrix
    inter_t = 31e-9  # Interphase thickness nm
    type_inter = 1  # 1- Soft interphase, 2- Hard interphase
    E_comp, nu_comp = computemechanicalprops(
        E_matrix, nu_matrix, E_cnt, nu_cnt, E_i, nu_i, L_cnt, D_cnt, inter_t, type_inter, vf
    )

    ###########################
    ###FRACTURE ENERGY###
    ###########################
    theta_min = 0  # [rad]
    theta_max = np.pi / 2  # [rad]
    p = 0.5
    q = 0.5
    L_min = 0
    L_max = L_cnt / 2  # Maximum Cnt length to integrate
    Gc = Energy_release_rate(
        theta_min,
        theta_max,
        L_min,
        L_max,
        vf,
        D_cnt,
        L_cnt,
        sigma_ult,
        Ac,
        tau_int,
        mu,
        E_cnt,
        G_m,
        p,
        q,
    )


    ###########################
    ###Electrical properties###
    ###########################
    dco = 1.870116073238336e-9  # Interparticle distance [nm]
    Lambdao = 0.500004882406769  # Height of the potential barrier [eV]
    # STRAIN SENSING CURVES
    # Maximum compression
    str_comp = -5
    # [%]
    # Maximum traction
    str_tens = 5
    # [%]
    # INPUT PARAMETERS
    # Matrix phase [Isotropic material]
    # densm = 1.12               # Density g/cm3
    sigma_m = 1.036000000000000e-10  # Electrical conductivity S/m
    # Electrical conductivities of CNTs to be analysed
    sigma_cnt = 10 ** 2
    # ELECTRICAL PROPERTIES
    # STRAIN SENSING CURVES
    # Maximum compression
    str_comp = -5
    # [%]
    # Maximum traction
    str_tens = 5
    # [%]
    sigmaEFF, L11_tract, L12_tract, L44_tract = electrical_prop(
        L_cnt, D_cnt, dco, Lambdao, sigma_m, sigma_cnt, vf, str_comp, str_tens
    )
    Eyoung = E_comp  # [Pa]
    Nu = nu_comp  # []
    Gc = Gc  # [J/m^2]
    cond = sigmaEFF[0, 0]  # [S/m]
    pi11 = L11_tract
    pi12 = L12_tract
    pi44 = L44_tract
    
    case_study(
        Vimp,
        Dimp,
        Inc,
        Eyoung,
        Nu,
        Gc,
        theta,
        cond,
        pi11,
        pi12,
        pi44,
        k_value,
        n_value,
        file_name,
        case_number,
        meshsize,
        meshsize2
        )
    post_process(file_name,file_name,case_number)
