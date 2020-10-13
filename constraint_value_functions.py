#Constraint Value Functions

import numpy as np


#%%PhI DP1 constraint equation
def phi_dp1(X,C,f_t):
    for i in X:
        if i.name == "body i":
            i_index=X.index(i)
        if i.name == "body j":
            j_index=X.index(i)            
            
    I=X[i_index]
    J=X[j_index]
    
    a_bar_i=C.a_bar_i
    a_bar_j=C.a_bar_j
    A_i=I.A_rotation
    A_j=J.A_rotation
    
    a_bar_i_T=np.transpose(a_bar_i)
    A_i_T=np.transpose(A_i)

    
    n1=np.dot(a_bar_i_T,A_i_T)
    n2=np.dot(n1,A_j)
    n3=np.dot(n2,a_bar_j)

    phi_dp1_value=n3-f_t

    return float(phi_dp1_value)
#%%
def phi_cd(X,C,ft):
    for i in X:
        if i.name == "body i":
            i_index=X.index(i)
        if i.name == "body j":
            j_index=X.index(i) 
               
    I=X[i_index]
    J=X[j_index]    
    

    A_i=I.A_rotation
    A_j=J.A_rotation

    r_j=J.q[0:3]
    s_bar_j=C.s_bar_j
    
    
    
    r_i=I.q[0:3]
    s_bar_i=C.s_bar_i
             
    c_T=np.transpose(C.c)
    
    if s_bar_j.all() == 0:
        d_ij=(-r_i-np.dot(A_i,s_bar_i))
    else:
        d_ij=(r_j+np.dot(A_j,s_bar_j)-r_i-np.dot(A_i,s_bar_i))
    
    PHI_CD=np.dot(c_T,d_ij)-ft
                        
    return float(PHI_CD)
#%%
def phi_dp2(X,C,ft):
    for i in X:
        if i.name == "body i":
            i_index=X.index(i)
        if i.name == "body j":
            j_index=X.index(i) 
               
    I=X[i_index]
    J=X[j_index]    
    
    a_bar_i=I.a_bar

    A_i=I.A_rotation
    A_j=J.A_rotation

    r_j=J.q[0:3]
    s_bar_j=J.s_bar
    
    r_i=I.q[0:3]
    s_bar_i=I.s_bar
    
    a_bar_i_T=np.transpose(a_bar_i)
    
    phi_dp2=np.dot(np.dot(a_bar_i_T,A_i),(r_j+np.dot(A_j,s_bar_j)-r_i-np.dot(A_i,s_bar_i)))
             

    return float(phi_dp2)
#%%
def phi_d(X,C,ft):
    for i in X:
        if i.name == "body i":
            i_index=X.index(i)
        if i.name == "body j":
            j_index=X.index(i) 
               
    I=X[i_index]
    J=X[j_index]    
    
    a_bar_i=I.a_bar

    A_i=I.A_rotation
    A_j=J.A_rotation

    r_j=J.q[0:3]
    s_bar_j=J.s_bar
    
    r_i=I.q[0:3]
    s_bar_i=I.s_bar
    

    a_bar_i_T=np.transpose(a_bar_i)
    d_ij=(r_j+np.dot(A_j,s_bar_j)-r_i-np.dot(A_i,s_bar_i))    
    
    phi_d=np.dot(np.transpose(d_ij),d_ij)-ft
             

    return float(phi_d)
