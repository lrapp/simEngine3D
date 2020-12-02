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
import xarray as xr
#%%



class sys():
    def __init__(self):    
        self.nb=0 #number of bodies
        self.bodies=[] #list of bodie objects
        self.constraints=[] #list of constraints
    

def dynamic_analysis(sys,t_start,t_step,t_end):
    SYS=sys()
    
    time_list=np.arange(t_start,t_end,t_step)
    
    SYS.bodies=data_file("C:\\Users\\Logan\\Desktop\\simEngine3D\\HW8_P1_revJoint.txt") #list of body objects
    SYS.constraints=constraints_in("C:\\Users\\Logan\\Desktop\\simEngine3D\\HW8_P1_revJoint.txt") #list of constraints

    SYS.nb=body_count(SYS.bodies)
    SYS.nc=len(SYS.constraints)
    

    
    
    #calculate intial position and velocity
    tol=1e-5
    check_phi(SYS.bodies,SYS.constraints,0,tol)    
    check_euler(SYS.bodies,SYS.constraints,0,tol)  
    
    compute_accel(SYS)
    
    outputs=build_DF(time_list) #construct DataFrame of length of time_list and columns defined in function

        
                                    
