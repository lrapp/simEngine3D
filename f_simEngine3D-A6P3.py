import os
os.chdir("C:\\Users\\Logan\\Desktop\\simEngine3D")
from pendulum_function import pendulum
import numpy as np
from constraints_in import constraints_in
from simEngine3D_dataload import data_file, DP1_PHI_partials,CD_PHI_partials,DP2_PHI_partials,D_PHI_partials, body
from simEngine3D_functions import build_p, build_A, calc_phi, calc_partials, build_ja,build_G
import matplotlib.pyplot as plt

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



F=np.array([0,0,-9.81])
F.shape=(3,1)

tau=np.array([0,0,0])
tau.shape=(3,1)

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
partial_list=[]
torque_list=[]
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
    n,partials,gamma_values=pendulum(X,constraint_list,time)
    partial_list.append(partials)
    body_list.append(n)
    x.append(X[1].q[0])
    y.append(X[1].q[1])
    z.append(X[1].q[2])
    r_dot.append(X[1].r_dot)
    r_d_dot.append(X[1].r_d_dot)    
    p_dot.append(X[1].p_dot)
    p_d_dot.append(X[1].p_d_dot)
    
    #Start inverse dynamics stuff
    G=build_G(X[1].q[3:])
    J_P=4*np.dot(np.dot(np.transpose(G),J_bar),G)    
    P=np.zeros([1,4])
    P[0,:]=np.transpose(X[1].q[3:])
    
    gamma_p=np.dot(P,X[1].p_d_dot)
    gamma_values=np.array(gamma_values)
    
    RHS=np.zeros([15,1])
    RHS[:3,0:1]=F
    RHS[7,0]=gamma_p
    RHS[8:,0]=gamma_values
     
    
    LHS=np.zeros([15,15])

    phi_r=partials[:,0:3]
    phi_p=partials[:,3:]
        
    LHS[0:3,0:3]=M
    LHS[0:3,8:]=np.transpose(phi_r)
    
    LHS[3:7,3:7]=J_P
    LHS[3:7,7:8]=np.transpose(P)
    LHS[3:7,8:]=np.transpose(phi_p)
    
    LHS[7,3:7]=P
       
    LHS[8:,0:3]=phi_r
    LHS[8:,3:7]=phi_p
       
    unknowns=np.dot(np.linalg.inv(LHS),RHS)
    
    rdd=unknowns[0:3]
    pdd=unknowns[3:7]
    lambda_P=unknowns[7]
    lagrange=unknowns[8:]
    
    torque=-np.dot(np.transpose(phi_p),lagrange)
    torque_list.append(torque)
    
        
    X=body_list[i]




    
#%%

#plt.plot()
#plt.axis([-1.5,1.5,-2,2])
#for i in range(0,len(y)):
#    plt.plot(float(y[i]),float(z[i]),'o')
#    plt.pause(0.001)
#%%

plt.plot(y,z)
#%%

torque_plot=[]
for i in torque_list:
    torque_plot.append(i[2])
    
plt.plot(tm,torque_plot)


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
