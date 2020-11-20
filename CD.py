# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 11:04:04 2020

@author: Logan
"""
import numpy as np

# calculate a_tilde
def tilde(a):
    ax, ay, az = a[0], a[1], a[2]
    t = np.array([[0, -az, ay],
                  [az, 0, -ax],
                  [-ay, ax, 0]])
    return t

#%%
def build_A(p):
        import numpy as np
        A=np.empty([3,3])
           
        e0=p[0]
        e1=p[1]
        e2=p[2]
        e3=p[3]
        
        A[0,0]=2*((e0**2+e1**2)-0.5)
        A[0,1]=2*(e1*e2-e0*e3)
        A[0,2]=2*(e1*e3+e0*e2)
        
        A[1,0]=2*(e1*e2+e0*e3)
        A[1,1]=2*((e0**2+e2**2)-0.5)        
        A[1,2]=2*(e2*e3-e0*e1)
        
        A[2,0]=2*(e1*e3-e0*e2)
        A[2,1]=2*(e2*e3+e0*e1)
        A[2,2]=2*((e0**2+e3**2)-0.5)        
        
        return A
#%%    
# get B matrix from p and a_bar
def getB(p, a_bar):
    B     = np.zeros((3, 4))
    e0, e = p[0], p[1:]
    # first column of B matrix
    c = e0 * np.identity(3) + tilde(e)
    # 3*3 matrix inside B
    m = np.matmul(c, tilde(a_bar))

    B[:, 0]  = np.dot(c, a_bar)
    B[:, 1:] = np.outer(e, a_bar) - m
    B *= 2.
    return B
#%%
#%%
def build_B(p_vector,a_vector):
    e0=p_vector[0]
    e=p_vector[1:]
    a_bar=a_vector
    I3=np.identity(3)
    
    a_bar_T=np.transpose(a_bar)
    a_bar_tilde=tilde(a_bar)
    e_tilde=tilde(e)
    
    B=np.zeros([3,4])
    X=(2*(np.dot(np.dot(float(e0),I3)+e_tilde,a_bar)))
       
    Y=(2*(np.dot(e,a_bar_T)-np.dot(np.dot(float(e0),I3)+e_tilde,a_bar_tilde)))
    for i in range(0,len(X)):
        B[i,0]=X[i]
    B[0,1:]=Y[0]
    B[1,1:]=Y[1]
    B[2,1:]=Y[2]

    return B

#%%
def CD_PHI_r(X,C):
    res = np.zeros(6)
    c=C.c
    c_T=np.transpose(c)
    j_ground=X[1].ground
    
    res[:3]=-c_T
    res[3:]=c_T
    if j_ground:
        return res[:3]
    else:
        return res
    #%%
    
def CD_PHI_p(X,C):
    res = np.zeros(8)
    c=C.c
    j_ground=X[1].ground
    
    s_bar_i=C.s_bar_i
    p=X[0].q[3:]

    res[:3]=-np.dot(c,build_B(X[0].q[3:],s_bar_i))
    res[3:]=c_T
    if j_ground:
        return res[:3]
    else:
        return res    
#%%
def CD_PHI_partials(X,C):
        #Set bodys i and j 
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
    

    p_i=I.q[3:]
    p_j=J.q[3:]
    
    s_bar_i=C.s_bar_i
    s_bar_j=C.s_bar_j
    

    
    p_i.shape=(4,1)
    s_bar_i.shape=(3,1)
    
    c=C.c
    c_T=np.transpose(c)
    
    first=-c_T
    second=np.dot(-c_T,build_B(p_i,s_bar_i))
    third=c_T
    fourth=np.dot(c_T,build_B(p_j,s_bar_j))
    
    if J.ground=='True':
        r_v = np.concatenate((first[0],second[0]))
    else:
        r_v= np.concatenate((first[0], third[0], second[0], fourth[0]))
        
    r_v.shape=(1,7)
    return r_v
    
#%%

def CD_phi(X,C,ft):
    body_1=C.body_i_name
    body_2=C.body_j_name
    
    for i in X:
        if i.name == body_1:
            i_index=X.index(i)
        if i.name == body_2:
            j_index=X.index(i) 
               
    I=X[i_index]
    J=X[j_index]    
    
    A_i=build_A(I.q[3:])
    A_j=build_A(J.q[3:])

    r_j=J.q[0:3]
    r_j.shape=(3)
    s_bar_j=C.s_bar_j
    
    r_i=I.q[0:3]
    r_i.shape=(3)
    s_bar_i=C.s_bar_i
             
    c_T=np.transpose(C.c)[0]
    
    s_j_Q=np.matmul(A_j,s_bar_j)
    s_j_Q.shape=(3)
    s_i_P=np.matmul(A_i,s_bar_i)
    s_i_P.shape=(3)

    PHI_CD=np.dot(c_T,r_j+s_j_Q-r_i-s_i_P)-ft
                        
    return float(PHI_CD)  

#%%
def CD_nue(X,C,dft):

    if C.name == "driving":
        dft=dft
    else:
        dft=0
        
    return dft

#%%
#%%
def gamma_CD(X,C,ddf):
        #Set bodys i and j 
    for i in X:
        if i.name == "body i":
            i_index=X.index(i)
        if i.name == "body j":
            j_index=X.index(i) 
            
    I=X[i_index]
    J=X[j_index]  
    
    p_dot_i=I.p_dot
    p_dot_j=J.p_dot
    
    p_dot_i.shape=(4,1)
    
    s_bar_i=C.s_bar_i
    s_bar_j=C.s_bar_j
    
    c=C.c
    
    c_T=np.transpose(c)

    if C.name == "driving":
        ddf=ddf
    else:
        ddf=0
        

    B_dot_i=build_B(p_dot_i,s_bar_i)
    
    if s_bar_j.all() == 0:
        B_dot_j=np.zeros([3,4])
    else:
        B_dot_j=build_B(p_dot_j,s_bar_j)
    
    gamma_CD=np.dot(np.dot(c_T,B_dot_i),p_dot_i)-np.dot(np.dot(c_T,B_dot_j),p_dot_j)+ddf
    
    return float(gamma_CD)