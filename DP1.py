# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 10:58:02 2020

@author: Logan
"""

from simEngine3D_functions import build_A,build_B
import numpy as np

#%%PhI DP1 constraint equation
def DP1_phi(X,C,f_t):
    body_1=C.body_i_name
    body_2=C.body_j_name
    
    for i in X:
        if i.name == body_1:
            i_index=X.index(i)
        if i.name == body_2:
            j_index=X.index(i) 
               
    I=X[i_index]
    J=X[j_index]    
    
    a_bar_i=C.a_bar_i
    a_bar_j=C.a_bar_j
    

    A_i=build_A(I.q[3:])
    A_j=build_A(J.q[3:])
    
    a_bar_i_T=np.transpose(a_bar_i)
    A_i_T=np.transpose(A_i)

    
    n1=np.dot(a_bar_i_T,A_i_T)
    n2=np.dot(n1,A_j)
    n3=np.dot(n2,a_bar_j)

    phi_dp1_value=n3-f_t

    return float(phi_dp1_value)

#%%
def DP1_PHI_partials(X,C):
    body_1=C.body_i_name
    body_2=C.body_j_name
    
    for i in X:
        if i.name == body_1:
            i_index=X.index(i)
        if i.name == body_2:
            j_index=X.index(i) 
               
    I=X[i_index]
    J=X[j_index]  
    #Set a and A for bodies i and j
    a_bar_i=C.a_bar_i
    a_bar_j=C.a_bar_j
    
    A_i=build_A(I.q[3:])
    A_j=build_A(J.q[3:])
    
    #Performe some transposes 
    a_bar_i_T=np.transpose(a_bar_i)
    a_bar_j_T=np.transpose(a_bar_j)
    A_i_T=np.transpose(A_i)
    A_j_T=np.transpose(A_j)

    p_i=I.q[3:]
    p_j=J.q[3:]
    p_dot_i=I.p_dot
    p_dot_j=J.p_dot

    a_j=np.dot(A_j,a_bar_j)
    a_i=np.dot(A_i,a_bar_i)
    a_i_T=np.transpose(a_i)


    a_j_T=np.transpose(a_j)
    
    
    p_i.shape=(4,1)
    a_bar_i.shape=(3,1)
    B_i=build_B(p_i,a_bar_i)
    B_j=build_B(p_j,a_bar_j)
    
    PHI_DP1_p_i=np.dot(a_j_T,B_i)
    PHI_DP1_p_j=np.dot(a_i_T,B_j)   
    
    PHI_DP1_p_i.shape=(1,4)
    PHI_DP1_p_j.shape=(1,4)
    
    PHI_DP1_r_i=np.array([0,0,0])
    PHI_DP1_r_j=np.array([0,0,0])
    
    PHI_DP1_r_i.shape=(1,3)
    PHI_DP1_r_j.shape=(1,3)    
    
    if J.ground=='True':
        return np.concatenate((PHI_DP1_r_i,PHI_DP1_p_i),axis=1)
    else:
        # return np.concatenate((PHI_DP1_r_i,PHI_DP1_r_j,PHI_DP1_p_i,PHI_DP1_p_j),axis=1)        
        return np.concatenate((PHI_DP1_r_i,PHI_DP1_p_i,PHI_DP1_r_j,PHI_DP1_p_j),axis=1)

#%%
def DP1_nue(X,C,dft):

    if C.name == "driving":
        dft=dft
    else:
        dft=0

    return dft
#%%    
def gamma_DP1(X,C,ddf):
        #Set bodys i and j 
    for i in X:
        if i.name == "body i":
            i_index=X.index(i)
        if i.name == "body j":
            j_index=X.index(i) 
            
    I=X[i_index]
    J=X[j_index]  
    
    #Set a and A for bodies i and j
    a_bar_i=C.a_bar_i
    a_bar_j=C.a_bar_j
    A_i=build_A(I.q[3:])
    A_j=build_A(J.q[3:])
    
    p_i=I.q[3:]
    p_j=J.q[3:]
    
    p_i.shape=(4,1)
    
    p_dot_i=I.p_dot
    p_dot_j=J.p_dot

    a_i=np.dot(A_i,a_bar_i)
    a_j=np.dot(A_j,a_bar_j)

    B_dot_i=build_B(p_dot_i,a_bar_i)    
    B_dot_j=build_B(p_dot_j,a_bar_j)
    
    B_i=build_B(p_i,a_bar_i)
    B_j=build_B(p_j,a_bar_j)
    
    a_dot_i=np.dot(B_i,p_dot_i)
    a_dot_j=np.dot(B_j,p_dot_j)
    
    a_i_T=np.transpose(a_i)
    a_j_T=np.transpose(a_j)
    
    a_dot_i_T=np.transpose(a_dot_i)
    
    if C.name == "driving":
        ddt=ddf
    else:
        ddt=0
    
    gamma_DP1=np.dot(np.dot(-a_i_T,B_dot_j),p_dot_j)-np.dot(np.dot(a_j_T,B_dot_i),p_dot_i)-2*np.dot(a_dot_i_T,a_dot_j)+ddt
    
    return float(gamma_DP1)
