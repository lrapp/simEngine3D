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

# get A matrix from euler parameters
def getA(p):
    e0, e1, e2, e3 = p[0], p[1], p[2], p[3]
    A = 2 * np.array([[e0**2 + e1**2 - 0.5, e1*e2 - e0*e3, e1*e3 + e0*e2],
                      [e1*e2 + e0*e3, e0**2 + e2**2 - 0.5, e2*e3 - e0*e1],
                      [e1*e3 - e0*e2, e2*e3 + e0*e1, e0**2 + e3**2 - 0.5]])
    return A

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
    A_i=I.A_rotation
    A_j=J.A_rotation
    
    #Performe some transposes 
    a_bar_i_T=np.transpose(a_bar_i)
    a_bar_j_T=np.transpose(a_bar_j)
    A_i_T=np.transpose(A_i)
    A_j_T=np.transpose(A_j)
    
    
    #define G and p_dot

    p_dot_i=I.p_dot
    p_dot_j=J.p_dot
    
    p_i=I.q[3:]
    p_j=J.q[3:]
    
    s_bar_i=C.s_bar_i
    s_bar_j=C.s_bar_j
    
    a_j=np.dot(A_j,a_bar_j)
    a_i=np.dot(A_i,a_bar_i)
    a_i_T=np.transpose(a_i)


    a_bar_tilde_i=tilde(a_bar_i)

    a_bar_tilde_j=tilde(a_bar_j)
    a_bar_tilde_j_T=np.transpose(a_bar_tilde_j)
    
    B_i=build_B(p_i,s_bar_i)
    B_j=build_B(p_j,s_bar_j)
    
    c=C.c
    c_T=np.transpose(c)
    
    first=-c_T
    second=np.dot(-c_T,B_i)
    third=c_T
    fourth=np.dot(c_T,B_j)
    
    
    return first[0], second[0], third[0], fourth[0]