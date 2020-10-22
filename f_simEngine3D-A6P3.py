import os
os.chdir("C:\\Users\\Logan\\Desktop\\simEngine3D")
from pendulum_function import pendulum
import numpy as np
from constraints_in import constraints_in
from simEngine3D_dataload import data_file, DP1_PHI_partials,CD_PHI_partials,DP2_PHI_partials,D_PHI_partials, body
from simEngine3D_functions import build_p, build_A, calc_phi, calc_partials, build_ja,build_G,tilde
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


F=np.array([0,0,m*-9.81])
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
    n,partials,gamma_values,pi_values=pendulum(X,constraint_list,time)
    partial_list.append(partials)
    body_list.append(n)
    x.append(X[1].q[0])
    y.append(X[1].q[1])
    z.append(X[1].q[2])
    r_dot.append(X[1].r_dot)
    r_d_dot.append(X[1].r_d_dot)    
    p_dot.append(X[1].p_dot)
    p_d_dot.append(X[1].p_d_dot)
    
    pi_values_a=np.zeros([6,3])
    for k in range(0,len(pi_values)):
        pi_values_a[k,:]=pi_values[k]
    
    #Start inverse dynamics stuff
    ## LHS is blue stuff on left of last equation on pg 478 of 2019 lecture notes
    ##RHS is blue stuff on right of last equation 
    
    phi_r=partials[:,0:3]
    phi_p=partials[:,3:]

    
    G=build_G(X[1].q[3:])
    pi_partials=0.5*np.dot(phi_p[0:6],G.T)    ##for r-omega
    
    J_P=4*np.dot(np.dot(np.transpose(G),J_bar),G)    
    P=np.zeros([1,4])
    P[0,:]=np.transpose(X[1].q[3:])
    
    gamma_p=np.dot(P,X[1].p_d_dot) #slide 475
#    gamma_p=-2*np.dot(X[1].p_dot.T,X[1].p_dot) #slide 473
    gamma_values=np.array(gamma_values)
    
    G_dot=build_G(X[1].p_dot)
#    tau_hat=8*np.dot(np.dot(np.dot(G_dot.T,J_bar),G_dot),X[1].q[3:])
    #Switching to r-omega formulation because I can't figure out the r-p form
    omega_bar=2*np.dot(G,X[1].p_dot)
    omega_bar_tilde=tilde(omega_bar)
    tau=np.dot(np.dot(omega_bar_tilde,J_bar),omega_bar)
    
    omega_bar_dot=2*np.dot(G,X[1].p_d_dot)

    LHS=np.zeros([6,6])
    LHS[0:3,0:6]=phi_r[0:6,:].T
    LHS[3:,0:6]=pi_partials[0:6,:].T
       
    RHS=np.zeros([6,1])
    RHS[0:3,0:1]=-(np.dot(M,X[1].r_d_dot)-F)
    RHS[3:6,0:1]=-(np.dot(J_bar,omega_bar_dot)-tau)
    
    unknowns=np.dot(np.linalg.inv(LHS),RHS)      
    
    #r-p formulation
#    LHS=np.zeros([8,8])
#    LHS[0:3,0:7]=phi_r.T    
#    LHS[3:7,0:7]=phi_p.T
#    LHS[3:7,7:8]=P.T
#    LHS[7:,4:]=P
#    
#    RHS=np.zeros([8,1])
#    RHS[:3,0:1]=F-np.dot(M,X[1].r_d_dot)
#    RHS[3:7]=tau_hat-np.dot(J_P,X[1].p_d_dot)
#    RHS[7]=gamma_p
#       
#    unknowns=np.dot(np.linalg.inv(LHS),RHS)    
#    
#    lagrange=unknowns[:7]
#    lambda_p=unknowns[7]
#    
#    torque=-np.dot(phi_p.T,lagrange)
#    torque_list.append(torque)

    torque=np.dot(-pi_partials[0:6,:].T,unknowns)
#    torque=np.dot(-pi_values[-1].T,unknowns[-1])
#    torque=np.dot(-pi_values_a.T,unknowns)
    torque_list.append(torque)
        
    X=body_list[i]




    
#%%

#plt.plot()
#plt.axis([-1.5,1.5,-2,2])
#for i in range(0,len(y)):
#    plt.plot(float(y[i]),float(z[i]),'o')
#    plt.pause(0.001)
#%% position of pendulum
plt.figure()
plt.plot(y,z)
#%%

torque_plot0=[]
torque_plot1=[]
torque_plot2=[]
#torque_plot3=[]

for i in torque_list:
    torque_plot0.append(i[0])
    torque_plot1.append(i[1])
    torque_plot2.append(i[2])
#    torque_plot3.append(i[3])
    
plt.figure()
#plt.plot(tm,torque_plot0,label='x)
#plt.plot(tm,torque_plot1)
plt.plot(tm,torque_plot2)
plt.title("Torque")
plt.ylabel("Torque (N*m)")
plt.xlabel("time")
#plt.plot(tm,torque_plot3)


#%%
plt.figure()
plt.plot(tm,x,label='x')
plt.plot(tm,y,label='y')
plt.plot(tm,z,label='z')
plt.title("x,y and z position vs time")
plt.legend()
plt.xlabel('time')

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
plt.title("Acceleration")
