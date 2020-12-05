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
from CD import *
from DP1 import *
from DP2 import *
from D import *
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#%%



class sys():
    def __init__(self):    
        self.nb=0 #number of bodies
        self.bodies=[] #list of bodie objects
        self.constraints=[] #list of constraints
    
    
def dynamic_analysis(SYS,t_start,t_step,t_end,file):

    SYS.h=t_step
    time_list=np.arange(t_start,t_end,t_step)
    
    SYS.bodies=data_file(file) #list of body objects
    SYS.constraints=constraints_in(file) #list of constraints

    SYS.nb=body_count(SYS.bodies)
    SYS.nc=len(SYS.constraints)
    
    SYS.outputs=build_DF(SYS,time_list) #construct DataFrame of length of time_list and columns defined in function
    
   #set time step index                  
    SYS.n=0
    SYS.time=time_list[SYS.n]

    #calculate intial position and velocity
    tol=1e-5
    check_phi(SYS.bodies,SYS.constraints,0,tol)    
    check_euler(SYS.bodies,SYS.constraints,0,tol) 
    
    #update time step with position and velocities
    update_level0(SYS)
    update_level1(SYS)

    compute_accel(SYS) #updates level2s (ddr and ddp)
    
    for i in range(1,len(time_list)):
        SYS.n=i
        BDF(SYS,SYS.n)

    return SYS

def inverse_dynamics_analysis(SYS,t_start,t_step,t_end,file):
    SYS.h=t_step
    num_steps=int((t_end-t_start) / t_step)
    time_list=np.arange(t_start,t_end,t_step)

    SYS.bodies=data_file(file) #list of body objects
    SYS.constraints=constraints_in(file) #list of constraints
    
    SYS.nb=body_count(SYS.bodies)
    SYS.nc=len(SYS.constraints)
    
    SYS.outputs=build_DF(SYS,time_list) #construct DataFrame of length of time_list and columns defined in function
    
    
    g=SYS.bodies[0].gravity
    m=SYS.bodies[0].mass
    
    SYS.torque=np.zeros((num_steps, 3))
    
    phi=calc_phi(SYS.bodies,SYS.constraints,0)
    phi_q=calc_partials(SYS.bodies,SYS.constraints)
    
    max_iter=20
    tol=1e-15
    step=0
    
    q_old=SYS.bodies[0].q.T[0]
    
    while step < num_steps:
         t=step*SYS.h
         # position analysis using Newton-Raphson
         diff = float("inf")
         n = 0
         SYS.n=n
         SYS.time=t
     
         while diff > tol and n < max_iter:
             q_new = q_old - np.matmul(np.linalg.inv(phi_q), phi)
             diff  = np.linalg.norm(q_new - q_old)
             q_old = q_new
             
             to_X=q_new.copy()
             to_X.shape=(7,1)
             SYS.bodies[0].q=to_X
         
             
             phi = calc_phi(SYS.bodies,SYS.constraints,t)
             phi_q = calc_partials(SYS.bodies,SYS.constraints)
     
             n += 1    
          
             
         r=to_X[:3]
         p=to_X[3:]
         
         nu = calc_nue(SYS)
         q_i_dot=np.linalg.solve(phi_q,nu)
        
         dr=q_i_dot[:3]
         dp=q_i_dot[3:]
                 
        
         gamma=calc_gamma(SYS)
         q_i_ddot = np.linalg.solve(phi_q,gamma)
        
         ddr=q_i_ddot[:3]
         ddp=q_i_ddot[3:]
         
         update_bodies(SYS,r,p,dr,dp,ddr,ddp) #update object with final values  
         update_level0(SYS) #update position in output object
         update_level1(SYS) #update velocity in output object
         update_level2(SYS) #update velocity in output object            
            
         #build matricies for inverse dynamics 
         M=build_M(SYS)
         F=build_F(SYS)
         J_p=build_JP(SYS)
         P=build_P(SYS)
         tau_hat=build_tau_hat(SYS)
                
         RHS_1=np.concatenate((F- np.matmul(M,ddr)))
         RHS_2=np.concatenate((tau_hat - np.matmul(J_p,ddp)))
        
         RHS=np.concatenate((RHS_1,RHS_2))
        
           
         lamb = np.linalg.solve(np.transpose(phi_q), RHS)
         Pi = 0.5 * np.matmul(phi_q[:,3:], np.transpose(getE(q_new[3:])))      
        
        
#          rhs = np.concatenate((g - np.matmul(m * np.eye(3), q_i_ddot[:3]), # r
#                               np.matmul(-getJp(q_new[3:], m, width, 2.* L), q_i_ddot[3:]))) # p
# #        
#          lamb = np.linalg.solve(np.transpose(phi_q), rhs)
#          Pi = 0.5 * np.matmul(phi_q[:,3:], np.transpose(getE(q_new[3:])))   
        
         SYS.torque[step,0] = -Pi[5,0] * lamb[5]
         SYS.torque[step,1] = -Pi[5,1] * lamb[3]
         SYS.torque[step,2] = -Pi[5,0] * lamb[4]
        
         step +=1        
    
    return SYS
    
    
