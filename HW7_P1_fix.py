import os
os.chdir("C:\\Users\\Logan\\Desktop\\simEngine3D")
from pendulum_function import pendulum
import numpy as np
from constraints_in import constraints_in
from simEngine3D_dataload import data_file, DP1_PHI_partials,CD_PHI_partials,DP2_PHI_partials,D_PHI_partials, body
from simEngine3D_functions import build_p, build_A, calc_phi, calc_partials, build_ja,build_G,tilde,build_p_from_A
import matplotlib.pyplot as plt

X=data_file("C:\\Users\\Logan\\Desktop\\simEngine3D\\revJoint_fix.txt")
constraint_list=constraints_in("C:\\Users\\Logan\\Desktop\\simEngine3D\\revJoint_fix.txt")

L=2
xa=0.05*0.05
volume=xa*L
rho=7800
m=rho*volume

M=m*np.identity(3)
J_bar=np.zeros([3,3])

b=0.05/2
c=0.05/2


J_bar[0,0]=1/12*m*(b**2+c**2)
J_bar[1,1]=1/12*m*(L**2+c**2)
J_bar[2,2]=1/12*m*(L**2+b**2)


F=np.array([0,0,m*-9.81]) 
F.shape=(3,1)

tau=np.array([0,0,0])
tau.shape=(3,1)


r_i=np.array([0, np.sin(np.pi/4.) * L, -np.cos(np.pi/4.) * L])
r_i.shape=(3,1)
X[0].q[0:3]=r_i

A = np.array([[0, 0, 1],
              [ np.sin(np.pi/4.), np.cos(np.pi/4.), 0],
              [-np.cos(np.pi/4.), np.sin(np.pi/4.), 0]])
p_i = build_p_from_A(A)
p_i_dot = np.array([0,0,0,0])
X[0].A_rotation=A
X[0].q[3:]=p_i 
X[0].ground=False
X[1].ground=True

phi=calc_phi(X,constraint_list,0)
phi_q=calc_partials(X,constraint_list)
