# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 10:54:23 2020

@author: Logan
"""
import numpy as np
from simEngine3D_functions import build_p, build_A, calc_phi, calc_partials, build_ja,build_G,tilde, build_B


def calc_partials_lagrange(X,cl,lagrange):
    phi_partial_lagrange=[]
    for i in cl:
            index=cl.index(i)
            if i.type.strip(" ' ") == 'CD':
                phi_lagrange=CD_phi_parital_lagrange(X,i,lagrange[index])
                phi_partial_lagrange.append(phi_lagrange)
            if i.type.strip(" ' ") == 'DP1':
                phi_lagrange=DP1_phi_parital_lagrange(X,i,lagrange[index])
                phi_partial_lagrange.append(phi_lagrange)    
    return phi_partial_lagrange


def K(a_bar,b):
    K=np.zeros([4,4])
    b_tilde=tilde(b)
    a_bar_tilde=tilde(a_bar)
    K[0,0]=a_bar.T @ b
    K[0,1:4] = a_bar.T @ b_tilde
    K[1:4,0:1]=a_bar_tilde @ b
    K[1:4,1:]=a_bar @ b.T + b @ a_bar.T - (a_bar.T @ b)*np.identity(3)
    return 2*K
     
def DP1_phi_parital_lagrange(X,C,lagrange):
    DP1_phi_partials=np.zeros([14,14])
    a_bar_i=C.a_bar_i
    a_bar_j=C.a_bar_j
    
    A_i=X[0].A_rotation
    A_j=X[1].A_rotation
    a_j=A_j @ a_bar_j
    a_i=A_i @ a_bar_i
    
    p_i=X[0].q[3:]
    p_j=X[1].q[3:]
    
    
    #assemble matrix from slide 634
    
    row3col4=build_B(p_i,a_bar_i).T @ build_B(p_j,a_bar_j)
    row4col3=build_B(p_j,a_bar_j).T @ build_B(p_i,a_bar_i)
    
    DP1_phi_partials[6:10,6:10]=K(a_bar_i,a_j)
    DP1_phi_partials[6:10,10:14]=row3col4
    DP1_phi_partials[10:14,6:10]=row4col3
    DP1_phi_partials[10:14,10:14]=K(a_bar_j,a_i)
    
    return lagrange*DP1_phi_partials
    
def CD_phi_parital_lagrange(X,C,lagrange):
    CD_phi_partials=np.zeros([14,14])
    s_bar_i=C.s_bar_i
    s_bar_j=C.s_bar_j
    c=C.c
    
    
    #assemble matrix from slide 634
    CD_phi_partials[6:10,6:10]=-K(s_bar_i,c)

    CD_phi_partials[10:14,10:14]=K(s_bar_j,c)
    
    return lagrange*CD_phi_partials
        