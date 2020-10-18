import os
os.chdir("C:\\Users\\Logan\\Desktop\\simEngine3D")
from pendulum_function import pendulum
import numpy as np
from constraints_in import constraints_in
from simEngine3D_dataload import data_file, DP1_PHI_partials,CD_PHI_partials,DP2_PHI_partials,D_PHI_partials, body
from simEngine3D_functions import build_p, build_A, calc_phi, calc_partials, build_ja
import matplotlib.pyplot as plt


X=data_file()
constraint_list=constraints_in()

p_j=np.array([0.7,0.1,0.7,0.1])
p_j.shape=(4,1)
X[1].A_rotation=build_A(p_j)

for i in range(0,len(X)):
    if X[i].name == "body j":
        X[i].q[3:]=p_j
           
X[1].q[:3,0]=np.array([0,1,-0.5])

t=np.arange(0,5,0.001)
body_list=[]
y=[]
z=[]
x=[]
for i in range(0,len(t)):
    time=t[i]
#    theta=np.pi/4*np.cos(2*time)
#    p_j=build_p(theta)
#    X[1].q[3:]=p_j
#    X[1].A_rotation=build_A(p_j)     
    X[1].A_rotation=build_A(X[1].q[3:])
    n=pendulum(X,constraint_list,time)
    body_list.append(n)
    x.append(X[1].q[0])
    y.append(X[1].q[1])
    z.append(X[1].q[2])
    X=body_list[i]

plt.plot(y,z)

#theta_list=[]
#for i in range(0,len(t)):
#    time=t[i]
#    theta_list.append(np.pi/4*np.cos(2*time))