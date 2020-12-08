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

tm=np.arange(0,10,0.01)
body_list=[]
y=[]
z=[]
x=[]
r_dot=[]
r_d_dot=[]
p_dot=[]
p_d_dot=[]
for i in range(0,len(tm)):
    time=tm[i]
 
    X[1].A_rotation=build_A(X[1].q[3:])
    n=pendulum(X,constraint_list,time)
    body_list.append(n)
    x.append(X[1].q[0])
    y.append(X[1].q[1])
    z.append(X[1].q[2])
    r_dot.append(X[1].r_dot)
    r_d_dot.append(X[1].r_d_dot)    
    p_dot.append(X[1].p_dot)
    p_d_dot.append(X[1].p_d_dot)

    X=body_list[i]
    


#plt.plot()
#plt.axis([-1.5,1.5,-2,2])
#for i in range(0,len(y)):
#    plt.plot(float(y[i]),float(z[i]),'o')
#    plt.pause(0.001)
#%%

plt.plot(y,z)


#theta_list=[]
#for i in range(0,len(t)):
#    time=t[i]
#    theta_list.append(np.pi/4*np.cos(2*time))

#%%
plt.figure()
plt.plot(tm,x)
plt.plot(tm,y)
plt.plot(tm,z)

#%%

r_dot_x=[]
r_dot_y=[]
r_dot_z=[]
for i in r_dot:
    r_dot_x.append(i[0])
    r_dot_y.append(i[1])
    r_dot_z.append(i[2])
    
plt.figure()
plt.title('Velocity of O')
plt.plot(tm,r_dot_x,label='r_dot_x')
plt.plot(tm,r_dot_y,label='r_dot_y')
plt.plot(tm,r_dot_z,label='r_dot_z')
plt.legend()

#%%

r_d_dot_x=[]
r_d_dot_y=[]
r_d_dot_z=[]
for i in r_d_dot:
    r_d_dot_x.append(i[0])
    r_d_dot_y.append(i[1])
    r_d_dot_z.append(i[2])
    
plt.figure()
plt.plot(tm,r_d_dot_x,label='r_d_dot_x')
plt.plot(tm,r_d_dot_y,label='r_d_dot_y')
plt.plot(tm,r_d_dot_z,label='r_d_dot_z')
