#simEngine3D-A6P2
import os
os.chdir("C:\\Users\\Logan\\Desktop\\simEngine3D")
from constraint_value_functions import phi_dp1,phi_dp2,phi_cd,phi_d
from nue_functions import DP1_nue,CD_nue,DP2_nue,D_nue
from gamma_functions import gamma_DP1,gamma_CD,gamma_DP2, gamma_D
from simEngine3D_dataload import data_file, DP1_PHI_partials,CD_PHI_partials,DP2_PHI_partials,D_PHI_partials, body
from HW6_Q2_flags import flags_in,points_in, f_in
from constraints_in import constraints_in
from simEngine3D_functions import build_p, build_A, calc_phi, calc_partials, build_ja,build_G
import numpy as np              
import sympy as sym


def pendulum(X,constraint_list,t):
#    f=np.cos(np.pi/4*np.cos(2*t))
#    df=np.cos(-1/2**np.pi*np.sin(2*t))
#    ddf=np.cos(-np.pi*np.cos(2*t))

    df=np.pi*np.sin(2*t)*np.cos((np.pi*np.cos(t*2))/4)*-1/2
    ddf=np.pi**2*np.sin(t*2)**2*np.sin((np.pi*np.cos(t*2))/4)*(-1/4)-np.pi*np.cos(t*2)*np.cos((np.pi*np.cos(t*2))/4)
                   
    
    phi_values=[]
    for i in constraint_list:
        if i.type.strip(" ' ") == 'CD':
            phi_values.append(phi_cd(X,i,eval(i.f)))
        elif i.type.strip(" ' ") == 'DP1':
            phi_values.append(phi_dp1(X,i,eval(i.f)))
    
    #Add euler parameter constraint
    phi_values.append(round(float((np.dot(np.transpose(X[1].q[3:]),X[1].q[3:])-1)),12))
    
    phi_partials_values=[]
    for i in constraint_list:
            if i.type.strip(" ' ") == 'CD':
                dri,dpi,drj,dpj=CD_PHI_partials(X,i)
                phi_partials_values.append([drj,dpj])
            if i.type.strip(" ' ") == 'DP1':
                dri,dpi,drj,dpj=DP1_PHI_partials(X,i)
                phi_partials_values.append([drj,dpj])

        
    #%% Assembel Jacobian
    ja_list=[]
    for i in range(0,len(phi_partials_values)):
        ca=np.zeros([7])
        #order for [d_ri d_pi d_rj d_pj]
        ca[0:3]=phi_partials_values[i][0]
        ca[3:7]=phi_partials_values[i][1]
    
        ca.shape=(1,7)
        ja_list.append(ca)
    
    
    jacobian=np.zeros([7,7])
    for i in range(0,len(ja_list)):
        jacobian[i,:]=ja_list[i]
        
    #add euler parameter normalization constraint to jacobian
    jacobian[6,3:]=2*X[1].q[3:].T
    #%%
    #Newton Raphson method
    q0=np.zeros([7,1])
    q0[:]=X[1].q
    
    phi=np.array(phi_values)    
    phi.shape=(7,1)
    
    error=10
    tol=1e-4
    counter=1
    
    while abs(error) > tol and counter < 30:
#        print(counter)
    #    q_new=q0-np.linalg.lstsq(jacobian,phi)[0]
        q_new=q0-np.dot(np.linalg.inv(jacobian),phi)
    #    q_new=q0-jacobian@phi
        error = np.linalg.norm(q0-q_new)
        
        X[1].q=q_new #Assing new q to body J
        X[1].A_rotation=build_A(q_new[3:])     #Build new rotation matrix from new p values
        phi=np.array(calc_phi(X,constraint_list,t)) #calculate new values of phi
        phi.shape=(7,1)
        partials_new=calc_partials(X,constraint_list) #calculate new partials
        jacobian=build_ja(X,partials_new) #build new jacobian
        q0=q_new
        counter=counter+1
    
    if counter == 30:
        print("Failed to converge")

                               
    nue_values=[]
    for i in constraint_list:
            if i.type.strip(" ' ") == 'CD':
                nue_values.append(CD_nue(X,i,df))
            if i.type.strip(" ' ") == 'DP1':
                nue_values.append(DP1_nue(X,i,df))            
                
    nue_values.append(0)
    
    #calculate p_dot
    q_dot=np.dot(np.linalg.inv(jacobian),nue_values)
    r_dot=q_dot[:3]
    r_dot.shape=(3,1)
    X[1].r_dot=r_dot
     
    p_dot=q_dot[3:]
    p_dot.shape=(4,1)
    X[1].p_dot=p_dot
     
    gamma_values=[]
    for i in constraint_list:
            if i.type.strip(" ' ") == 'CD':
                gamma_values.append(gamma_CD(X,i,ddf))
            if i.type.strip(" ' ") == 'DP1':
                gamma_values.append(gamma_DP1(X,i,ddf))
    
        
    gamma_values.append(float(-2*np.dot(np.transpose(X[1].p_dot),X[1].p_dot)))
    
    q_d_dot=np.dot(np.linalg.inv(jacobian),gamma_values)
    r_d_dot=q_d_dot[:3]
    r_d_dot.shape=(3,1)
    X[1].r_d_dot=r_d_dot
    p_d_dot=q_d_dot[3:]
    p_d_dot.shape=(4,1)
    X[1].p_d_dot=p_d_dot
   
    return X,jacobian,gamma_values
    