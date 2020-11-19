# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 14:27:42 2020
Restarting with pieces of my SimEngine3D from homework.  Starting over here for a clean slate.

@author: Logan
"""
import os
os.chdir("C:\\Users\\Logan\\Desktop\\simEngine3D")
from simEngine3D_dataload import data_file, constraints_in
from simEngine3D_functions import *


import numpy as np

file="C:\\Users\\Logan\\Desktop\\simEngine3D\\revJoint_fix2.txt"
#%%
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


    r_i = np.array([0, np.sin(np.pi/4.) * L, -np.cos(np.pi/4.) * L])
    r_i.shape=(3,1)

    X[0].q[:3]=r_i
    A = np.array([[0, 0, 1],
                  [ np.sin(np.pi/4.), np.cos(np.pi/4.), 0],
                  [-np.cos(np.pi/4.), np.sin(np.pi/4.), 0]])
    
    X[0].q[3:]=build_p_from_A(A)
     

     
    phi=calc_phi(X,constraint_list,0)
    phi_q=calc_partials(X,constraint_list)

    # inverse dynamics settings
    t_start, t_end = 0, 10
    dt = 1e-3
    num_steps = int((t_end-t_start) / dt)
    step = 0
    max_iter = 10
    tol = 1e-15
    
    # create arrays to save time evolution results
    # torque
    torque = np.zeros((num_steps, 3))

    q_old = X[0].q.T[0]

    while step < num_steps:
        t=step*dt
        # position analysis using Newton-Raphson
        diff = float("inf")
        iter = 0

        while diff > tol and iter < max_iter:
            q_new = q_old - np.matmul(np.linalg.inv(phi_q), phi)
            diff  = np.linalg.norm(q_new - q_old)
            q_old = q_new
            
            phi = calc_phi(X,constraint_list,t)
            phi_q = calc_partials(X,constraint_list)


            iter += 1        
            
        nu = calc_nue(X,constraint_list,t)
        q_i_dot=np.linalg.solve(phi_q,nu)
        
        X[0].q=q_new
        X[0].r_dot=q_i_dot[:3]
        X[0].p_dot=q_i_dot[3:]

        
        gamma=calc_gamma(X,constraint_list,t)
