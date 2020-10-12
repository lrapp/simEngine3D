#simEngine3D-A6P2
import os
os.chdir("C:\\Users\\Logan\\Desktop\\ME-751\\HW6\\Q1")
from constraint_value_functions import phi_dp1,phi_dp2,phi_cd,phi_d
from nue_functions import DP1_nue,CD_nue,DP2_nue,D_nue
from gamma_functions import gamma_DP1,gamma_CD,gamma_DP2, gamma_D
from HW6_Q2_dataload import data_file, DP1_PHI_partials,CD_PHI_partials,DP2_PHI_partials,D_PHI_partials, body
from HW6_Q2_flags import flags_in, constraints_in,points_in, f_in

import numpy as np
J=body()
I=body()
#file="C:\\Users\\Logan\\Desktop\\ME-751\\HW6\\input.txt" 

theta=0
L=2
t=0

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



I.name="body i"
J.name="body j"
I.s_bar=np.array([-2,0,0])
I.q[:3]=np.array([0,L*np.cos(theta),-L*np.sin(theta)])
I.q[3:]=np.array([e0,e1,e2,e3])
ft=np.pi/4*np.cos(2*t)

c=[1,0,0]

X=[I,J]
cd1=phi_cd(X,ft,np.array([1,0,0]))
cd2=phi_cd(X,ft,np.array([0,1,0]))
cd3=phi_cd(X,ft,np.array([0,0,1]))

perp_1_bodies=[body(),body()]
perp_1_bodies[0].name="body i"
perp_1_bodies[1].name="body j"

perp_1_bodies[0].q[:3]=np.array([0,L*np.cos(theta),-L*np.sin(theta)])
perp_1_bodies[0].q[3:]=np.array([e0,e1,e2,e3])
perp_1_bodies[1].q[:3]=np.array([0,0,0])
perp_1_bodies[0].a_bar=np.array([1,0,0])

perp_1_bodies[1].a_bar=np.array([1,0,0])

dp1_1=phi_dp1(perp_1_bodies,ft)             