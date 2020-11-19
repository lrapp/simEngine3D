# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 14:27:42 2020
Restarting with pieces of my SimEngine3D from homework.  Starting over here for a clean slate.

@author: Logan
"""
import os
os.chdir("C:\\Users\\Logan\\Desktop\\simEngine3D")
from simEngine3D_dataload import data_file, constraints_in
from simEngine3D_functions import get_p, calc_phi

from CD import CD_phi
import numpy as np

file="C:\\Users\\Logan\\Desktop\\simEngine3D\\revJoint_fix2.txt"

def main(file):
    t = 0.
    L = 2.
    width = 0.05
    rho   = 7800.
    m     = 2.* L * width**2 * rho
    g     = np.array([0, 0, -m * 9.81])
    
    X=data_file(file) #list of body objects
    constraint_list=constraints_in(file) #list of constraints
    
    j_ground=True
    
    #setup G-RF                            
    x=np.array([1,0,0])
    y=np.array([0,1,0])
    z=np.array([0,0,1])
    
    #setup L-RF
    x_p = np.array([1,0,0])
    y_p = np.array([0,1,0])
    z_p = np.array([0,0,1])

    # settings for CD constraints
    # position of Q in G-RF, and it's fixed
    s_j_Q_bar = np.array([0,0,0])
    r_j = np.array([0,0,0])
    p_j = np.array([1,0,0,0])
    r_j_dot = np.array([0,0,0])
    p_j_dot = np.array([0,0,0,0])

    # position of Q in L-RF, consider theta = pi/4 as initial position
    # distance between O' and Q are fixed
    s_i_P_bar = np.array([-L,0,0])
    r_i = np.array([0, np.sin(np.pi/4.) * L, -np.cos(np.pi/4.) * L])
    r_i.shape=(3,1)
    r_i_dot = np.array([0,0,0])
    A = np.array([[0, 0, 1],
                  [ np.sin(np.pi/4.), np.cos(np.pi/4.), 0],
                  [-np.cos(np.pi/4.), np.sin(np.pi/4.), 0]])
    p_i = get_p(A)
    p_i.shape=(4,1)
    p_i_dot = np.array([0,0,0,0])
    
    X[0].q[:3]=r_i
    X[0].q[3:]=p_i
     
    phi=calc_phi(X,constraint_list,0)
    partials=calc_partials(X,constraint_list)
