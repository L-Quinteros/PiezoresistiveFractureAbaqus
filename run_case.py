from Function_study_case import *

# Parameters
###################################################################################
###################################################################################
# Volume fraction to study
vf = 0.01  # 1 % volume fraction
theta = 0.0  # Crack angle #case 1 and 2
case_number = 1
Vimp = 10.0  # [v] Applied Voltage
Dimp = 0.0001  # [m] Applied Displacement
Inc = 0.005  # Increment of steps
k_value = 50  #
n_value = 6  #
meshsize1_serie = [0.001]
meshsize2_serie = [1,2,3,4,5]


for idx1,meshsize1 in enumerate(meshsize1_serie):
    for idx2,meshsize2 in enumerate(meshsize2_serie):
        length_l = meshsize2*meshsize1
        file_name = 'case_length_scale'+length_l+'mesh_size'+meshsize1+'theta_'+theta+'k_'+k_value+'n_'+n_value
        run_the_case(vf,theta,Inc,case_number,Vimp,Dimp,file_name,k_value,n_value,meshsize1,length_l)
