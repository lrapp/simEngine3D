# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 20:17:53 2020

@author: Logan
"""
import numpy as np
from simEngine3D_functions import build_G

L=2
xa=0.05*0.05
volume=xa*L
rho=7800
m=rho*volume

M=m*np.identity(3)
J=np.zeros([3,3])

b=0.05/2
c=0.05/2


J[0,0]=1/12*m*(b**2+c**2)
J[1,1]=1/12*m*(L**2+c**2)
J[2,2]=1/12*m*(L**2+b**2)

