# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 16:59:00 2020

@author: Logan
"""

import numpy as np

theta=np.pi/4
A=np.zeros([3,3])
A[0,2]=-1
A[1,0]=np.sin(theta)
A[1,1]=-np.sin(theta)
A[2,0]=-np.cos(theta)
A[2,1]=-np.sin(theta)

s_bar=np.array([1,0,0])
s=np.dot(A,s_bar)

e0=np.sqrt((np.trace(A)+1)/4)
e1=(A[2,1]-A[1,2])/(4*e0)
e2=(A[0,2]-A[2,0])/(4*e0)
e3=(A[1,0]-A[0,1])/(4*e0)
