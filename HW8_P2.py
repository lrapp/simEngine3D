import os
os.chdir("C:\\Users\\Logan\\Desktop\\simEngine3D")
from pendulum_function import pendulum
import numpy as np
from constraints_in import constraints_in
from simEngine3D_dataload import data_file, DP1_PHI_partials,CD_PHI_partials,DP2_PHI_partials,D_PHI_partials, body
from simEngine3D_functions import build_p, build_A, calc_phi, calc_partials, build_ja,build_G,tilde,calc_nue,check_phi,calc_gamma,build_E,build_p_from_A,calc_partials_HW82
from simEngine3D_functions import calc_phi_HW82
import matplotlib.pyplot as plt
from new_partials import DP1_phi_parital_lagrange,CD_phi_parital_lagrange
import time as ttime

import time as ttime

tic=ttime.perf_counter()

# Read in body and constraint information
file="C:\\Users\\Logan\\Desktop\\simEngine3D\\HW8_P2_input.txt"
X=data_file(file)
constraint_list=constraints_in(file)

#
h=0.01
L=2
xa=0.05*0.05
volume=xa*L
rho=7800
m1=rho*volume
m2=rho*volume/2

counter=0
time_list=np.arange(0,200+h,h)
                       
time=time_list[0]

#df=np.pi*np.sin(2*time)*np.cos((np.pi*np.cos(time*2))/4)*-1/2
#ddf=np.pi**2*np.sin(time*2)**2*np.sin((np.pi*np.cos(time*2))/4)*(-1/4)-np.pi*np.cos(time*2)*np.cos((np.pi*np.cos(time*2))/4)
df=0
ddf=0

nb=2
nc=len(constraint_list)

r_list=np.zeros([3*nb,len(time_list)])
r_d_list=np.zeros([3*nb,len(time_list)])
r_dd_list=[]

p_list=np.zeros([4*nb,len(time_list)])
p_d_list=np.zeros([4*nb,len(time_list)])
p_dd_list=[]


lambda_p_list=[]
lagrange_list=[]
         
F1=np.array([0,0,m1*-9.81]) 
F1.shape=(3,1)
F2=np.array([0,0,m2*-9.81]) 
F2.shape=(3,1)

F=np.zeros([6,1])
F[0:3]=F1
F[3:6]=F2

tau=np.array([0,0,0])
tau.shape=(3,1)

# set body j starting values
r_start=np.array([0,L,0]) 
A_start=np.zeros([3,3])
A_start[0,2]=np.cos(0)
A_start[1,0]=np.cos(0)
A_start[2,1]=np.cos(0)
p_start=build_p_from_A(A_start)

r_start.shape=(3,1)
p_start.shape=(4,1)
X[1].q[0:3]=r_start
X[1].q[3:]=p_start 
X[1].A_rotation=build_A(X[1].q[3:])
X[1].r_dot=np.array([0,0,0])
X[1].p_dot=np.array([0,0,0,0])

r_start=np.array([0,2*L,-L/2])
A_start=np.zeros([3,3])
A_start[0,2]=np.cos(0)
A_start[1,1]=np.cos(0)
A_start[2,0]=np.cos(np.pi)
p_start=build_p_from_A(A_start)
r_start.shape=(3,1)
p_start.shape=(4,1)
X[2].q[0:3]=r_start 
X[2].q[3:]=p_start 
X[2].A_rotation=build_A(X[2].q[3:])
X[2].r_dot=np.array([0,0,0])
X[2].p_dot=np.array([0,0,0,0])

#
tol=1e-5
phi=calc_phi_HW82(X,constraint_list,0)
check_phi(X,constraint_list,0,tol)
PHI_P=np.zeros([2,1])
for i in range(1,len(X)):
    PHI_P[i-1]=1/2*(X[i].q[3:].T @ X[i].q[3:]) - 1/2
         
if PHI_P.any() > tol:
    print("initial conditions do not satisfy PHI_P=0")
else:
    print("initial conditions satisfy PHI_P=0")
    
    
partials=calc_partials_HW82(X,constraint_list)
jacobian=np.vstack(partials)

phi_r=np.zeros([nc,6])
phi_r[:5,:]=jacobian[0:5,0:6]
phi_r[5:,:]=jacobian[5:,3:9]

#phi_p=jacobian[:,13:]
phi_p=np.zeros([nc,8])
phi_p[:5,:]=jacobian[0:5,9:17]
phi_p[5:,:]=jacobian[5:,13:]

r_list[0:3,0:1]=X[1].q[0:3]
r_list[3:6,0:1]=X[2].q[0:3]

p_list[0:4,0:1]=X[1].q[3:]
p_list[4:8,0:1]=X[2].q[3:]

for i in X:
    i.r_dot.shape=(3,1)
    i.p_dot.shape=(4,1)
    
r_d_list[0:3,0:1]=X[1].r_dot
r_d_list[3:6,0:1]=X[2].r_dot
        
        
p_d_list[0:4,0:1]=X[1].p_dot
p_d_list[4:10,0:1]=X[2].p_dot

nue=calc_nue(X,constraint_list,0)

check_nue=phi_r @ r_d_list[:,0] + phi_p @ p_d_list[:,0]

if (check_nue == nue).any():
    print("initial conditions satisfy nue=0")
else:
    print("initial conditions do not satisfy nue=0") 
    
P=np.zeros([2,8])
P[0,0:4]=X[1].q[3:].T
P[1,4:8]=X[2].q[3:].T
 
check_P=P @ p_d_list[:,0]

if (check_P == np.array([0,0])).any():
    print("initial conditions satisfy P_p_dot=0")
else:
    print("initial conditions do not satisfy P_p_dot=0") 

r_dd_list=[]
p_dd_list=[]
lambda_p_list=[]
lagrange_list=[]
nue_list=[]
#%%

partials=calc_partials_HW82(X,constraint_list)
jacobian=np.vstack(partials)

phi_r=np.zeros([nc,6])
phi_r[:5,:]=jacobian[0:5,0:6]
phi_r[5:,:]=jacobian[5:,3:9]

#phi_p=jacobian[:,13:]
phi_p=np.zeros([nc,8])
phi_p[:5,:]=jacobian[0:5,9:17]
phi_p[5:,:]=jacobian[5:,13:]

nue_values=calc_nue(X,constraint_list,df)

#phi_r=np.zeros([10,3*nb])
#phi_r[:,0:3]=jacobian[:,0:3] #get ri
#phi_r[:,3:6]=jacobian[:,3:6]
#
#phi_p=np.zeros([10,8])
#phi_p[:,0:4]=jacobian[:,6:10]
#phi_p[:,4:8]=jacobian[:,10:14]

gamma_values=calc_gamma(X,constraint_list,ddf)
           
           
M1=m1*np.identity(3)
M2=m2*np.identity(3)

M=np.zeros([6,6])
M[0:3,0:3]=M1
M[3:6,3:6]=M2
 
J_bar1=np.zeros([3,3])

b=0.05/2
c=0.05/2
J_bar1[0,0]=1/12*m1*(b**2+c**2)
J_bar1[1,1]=1/12*m1*(L**2+c**2)
J_bar1[2,2]=1/12*m1*(L**2+b**2)

J_bar2=np.zeros([3,3])
J_bar2[0,0]=1/12*m1*(b**2+c**2)
J_bar2[1,1]=1/12*m1*(L/2**2+c**2)
J_bar2[2,2]=1/12*m1*(L/2**2+b**2)

J_bar=[J_bar1,J_bar2]

G=[]
for i in range(1,len(X)):
    G.append(build_G(X[i].q[3:]))

J_P_list=[]
for i in range(0,len(G)):
    J_P_list.append(4*np.dot(np.dot(np.transpose(G[i]),J_bar[i]),G[i]))    
    
J_P=np.zeros([8,8])
J_P[0:4,0:4]=J_P_list[0]
J_P[4:8,4:8]=J_P_list[1]

G_dot=[]
for i in range(1,len(X)):
    G_dot.append(build_G(X[i].p_dot))

tau_hat_list=[]
for i in range(0,len(G_dot)):
    tau_hat_list.append(8*np.dot(np.dot(np.dot(G_dot[i].T,J_bar[i]),G_dot[i]),X[i+1].q[3:]))

tau_hat=np.zeros([8,1])
tau_hat[0:4]=tau_hat_list[0]
tau_hat[4:8]=tau_hat_list[1]

P=np.zeros([2,8])
P[0,0:4]=X[1].q[3:].T
P[1,4:8]=X[2].q[3:].T

LHS=np.zeros([26,26])
LHS[0:6,0:6]=M
LHS[0:6,16:26]=phi_r.T 
   
LHS[6:14,6:14]=J_P
LHS[6:14,14:16]=P.T
LHS[6:14,16:26]=phi_p.T

LHS[14:16,6:14]=P
   
LHS[16:26,0:6]=phi_r
LHS[16:26,6:14]=phi_p
   
RHS=np.zeros([26,1])

gamma_p=[]
for i in range(1,len(X)):   
    gamma_p.append(-2*np.dot(X[i].p_dot.T,X[i].p_dot)) #slide 473

RHS[0:6]=F
RHS[6:14]=tau_hat
RHS[14]=gamma_p[0]
RHS[15]=gamma_p[1]
gamma_hat=np.array(gamma_values) #remove eulear parameterization constraint
gamma_hat.shape=(10,1)
RHS[16:26]=gamma_hat
   
unknowns=np.dot(np.linalg.inv(LHS),RHS) 
r_dd_list.append(unknowns[0:6])
p_dd_list.append(unknowns[6:14])
lambda_p_list.append(unknowns[14:16])
lagrange_list.append(unknowns[16:26])
  




#%%
n=1
time=time_list[n]
for ii in range(1,len(time_list)):
    
    time=time_list[ii]
    r_dd=r_dd_list[n-1]
    p_dd=p_dd_list[n-1]
    lagrange=lagrange_list[n-1]
    lambda_p=lambda_p_list[n-1]
    
        
    z_0=np.zeros([7*nb+1*nb+nc,1])
    z_0[0:3*nb]=r_dd_list[n-1]
    z_0[3*nb:7*nb]=p_dd_list[n-1]
    z_0[7*nb:7*nb+nb]=lambda_p_list[n-1]
    z_0[7*nb+nb:]=lagrange_list[n-1]
    
    tol=1e-4
    error=1
    count=1
    while abs(error) > tol and count < 30:
            if count == 1:
                z=z_0
            else:
                z=z
                
                
            beta_0=1
            alpha1=1
        
            lambda_p=z[14:16]             
            lagrange=z[16:26]
            
            C_r_n=alpha1*r_list[:,n-1]
            C_r_n.shape=(3*nb,1)
            C_p_n=alpha1*p_list[:,n-1]
            C_p_n.shape=(4*nb,1)
            
            C_r_dot_n=alpha1*r_d_list[:,n-1]
            C_r_dot_n.shape=(3*nb,1)
            C_p_dot_n=alpha1*p_d_list[:,n-1]
            C_p_dot_n.shape=(4*nb,1)         
        
                  
            r=C_r_n+(beta_0**2)*(h**2)*z[0:3*nb]
            p=C_p_n+(beta_0**2)*(h**2)*z[3*nb:7*nb]        
        
            r_d=C_r_dot_n+beta_0*h*z[0:3*nb]
            p_d=C_p_dot_n+beta_0*h*z[3*nb:7*nb]          
            
            
            X_n=X
            X_n[1].q[0:3]=r[0:3]
            X_n[2].q[0:3]=r[3:6]    
            X_n[1].q[3:]=p[0:4]
            X_n[2].q[3:]=p[4:8]
            X_n[1].A_rotation=build_A(X_n[1].q[3:])
            X_n[2].A_rotation=build_A(X_n[2].q[3:])
            
            X_n[1].r_dot=r_d[0:3]
            X_n[2].r_dot=r_d[3:6]
            
            X_n[1].p_dot=p_d[0:4]
            X_n[2].p_dot=p_d[4:8]    
            
            X_n[1].r_d_dot=z[0:3]
            X_n[2].r_d_dot=z[3:6]
            
            X_n[1].p_d_dot=z[6:10]
            X_n[2].p_d_dot=z[10:14]
            
            PHI=calc_phi_HW82(X_n,constraint_list,0)
            PHI=np.array(PHI)
            
            partials=calc_partials_HW82(X_n,constraint_list)
            jacobian=np.vstack(partials)
            #jacobian=build_ja(X,partials)
            
            phi_r=jacobian[:,3:9]
            
            phi_p=jacobian[:,13:]    
            
            nue_values=calc_nue(X_n,constraint_list,df)
            
            #jacobian=jacobian[:-2,:] #remove euler parameterization constraint
            
        #    phi_r=np.zeros([10,6])
        #    phi_r[:,0:3]=jacobian[:,0:3] #get ri
        #    phi_r[:,3:6]=jacobian[:,3:6]
        #    
        #    phi_p=np.zeros([10,8])
        #    phi_p[:,0:4]=jacobian[:,6:10]
        #    phi_p[:,4:8]=jacobian[:,10:14]    
        #    
            for i in range(1,len(X)):
                PHI_P[i-1]=1/2*(X_n[i].q[3:].T @ X_n[i].q[3:]) - 1/2
                     
            G_dot=[]
            for i in range(1,len(X)):
                G_dot.append(build_G(X_n[i].p_dot))
            
            tau_hat_list=[]
            for i in range(0,len(G_dot)):
                tau_hat_list.append(8*np.dot(np.dot(np.dot(G_dot[i].T,J_bar[i]),G_dot[i]),X_n[i+1].q[3:]))  
                    
            G=[]
            for i in range(1,len(X)):
                G.append(build_G(X_n[i].q[3:]))
            
            J_P_list=[]
            for i in range(0,len(G)):
                J_P_list.append(4*np.dot(np.dot(np.transpose(G[i]),J_bar[i]),G[i]))            
                
            P=np.zeros([2,8])
            P[0,0:4]=X[1].q[3:].T
            P[1,4:8]=X[2].q[3:].T
                    
             
            g=np.zeros([8*nb+nc,1])
            g[0:3*nb]=M @ z[0:3*nb] + phi_r.T @ z[16:26] - F
            g[3*nb:3*nb+4*nb] = J_P @ z[3*nb:3*nb+4*nb] + phi_p.T @ z[16:26] + P.T @ z[14:16]- tau_hat
            g[3*nb+4*nb:3*nb+4*nb+nb]=1/((beta_0**2)*(h**2))*PHI_P
            last_g=(1/((beta_0**2)*(h**2)) * PHI)
            last_g.shape=(nc,1)
            g[3*nb+4*nb+nb:3*nb+4*nb+nb+nc]=last_g
             
            
            PSI=np.zeros([8*nb+nc,8*nb+nc])
            PSI[0:3*nb,0:3*nb]=M
            PSI[0:3*nb,3*nb+4*nb+nb:3*nb+4*nb+nb+nc]=phi_r.T
            
            PSI[3*nb:3*nb+4*nb,3*nb:3*nb+4*nb]=J_P
            PSI[3*nb:3*nb+4*nb,3*nb+4*nb:3*nb+4*nb+nb]=P.T
            PSI[3*nb:3*nb+4*nb,3*nb+4*nb+nb:3*nb+4*nb+nb+nc]=phi_p.T
        
            PSI[3*nb+4*nb:3*nb+4*nb+nb,3*nb:3*nb+4*nb]=P
               
            PSI[3*nb+4*nb+nb:3*nb+4*nb+nb+nc,0:3*nb]=phi_r
            PSI[3*nb+4*nb+nb:3*nb+4*nb+nb+nc,3*nb:3*nb+4*nb]=phi_p
            
            delta_z=-np.dot(np.linalg.inv(PSI),g)           
            z=z+delta_z
            error=np.linalg.norm(delta_z)
            count=count+1
            if count > 29:
                print("did not converge")    
   
    nue_values=calc_nue(X_n,constraint_list,df)
    nue_list.append(nue_values)
    vel=phi_r @ r_d + phi_p @ p_d
#    vel_violation.append(np.linalg.norm(vel))
    r_list[:,n:n+1]=r
    p_list[:,n:n+1]=p

    r_d_list[:,n:n+1]=r_d
    p_d_list[:,n:n+1]=p_d
            
    r_dd_list.append(z[0:6])
    p_dd_list.append(z[6:14])
    lambda_p_list.append(z[14:16])
    lagrange_list.append(z[16:26])    

#    G=build_G(p)
#    torque=[]
#    for kk in range(0,len(constraint_list)):
#        torque.append(-1/2*G @ phi_p[kk,:] * z[8+kk])
        
#    torque=-0.5*(G @ phi_p[5] * z[13:14])
#    torque=phi_p.T @ z[8:14] + P.T @ z[7:8]
#    torque_list.append(torque)
    n=n+1
    
toc=ttime.perf_counter()
elapsed_time=toc-tic
print(elapsed_time)    
#%%

xx=r_list[0,:]
yy=r_list[1,:]
zz=r_list[2,:]
fig, axs = plt.subplots(3)
axs[0].set(title="Body 1 point O'")
axs[0].plot(time_list,xx)
axs[0].set(ylabel='x')
axs[1].plot(time_list,yy)
axs[1].set(ylabel='y')
axs[2].plot(time_list,zz)
axs[2].set(ylabel='z')
axs[2].set(xlabel='time')
#%%calculate omega
E_list=[]
for i in range(0,max(p_list.shape)):
    E_list.append(build_E(p_list[:4,i]))

omega_list=np.zeros([3,max(p_list.shape)])
for i in range(0,max(p_d_list.shape)):
    omega=2*(np.dot(E_list[i], p_d_list[:4,i]))
    omega.shape=(3,1)
    omega_list[:,i:i+1]=omega

omx=omega_list[0,:]
omy=omega_list[1,:]
omz=omega_list[2,:]
fig, axs = plt.subplots(3)
axs[0].set(title="Body 1 point O'")
axs[0].plot(time_list,omx)
axs[0].set(ylabel='omega x')
axs[1].plot(time_list,omy)
axs[1].set(ylabel='omega y')
axs[2].plot(time_list,omz)
axs[2].set(ylabel='omega z')
axs[2].set(xlabel='time')
#%%

xx2=r_list[3,:]
yy2=r_list[4,:]
zz2=r_list[5,:]
fig, axs = plt.subplots(3)
axs[0].set(title="Body 2 point O'")
axs[0].plot(time_list,xx2)
axs[0].set(ylabel='x')
axs[1].plot(time_list,yy2)
axs[1].set(ylabel='y')
axs[2].plot(time_list,zz2)
axs[2].set(ylabel='z')
axs[2].set(xlabel='time')
#%%calculate omega
E_list=[]
for i in range(0,max(p_list.shape)):
    E_list.append(build_E(p_list[4:8,i]))

omega_list=np.zeros([3,max(p_list.shape)])
for i in range(0,max(p_d_list.shape)):
    omega=2*(np.dot(E_list[i], p_d_list[4:8,i]))
    omega.shape=(3,1)
    omega_list[:,i:i+1]=omega

omx=omega_list[0,:]
omy=omega_list[1,:]
omz=omega_list[2,:]
fig, axs = plt.subplots(3)
axs[0].set(title="Body 2 point O'")
axs[0].plot(time_list,omx)
axs[0].set(ylabel='omega x')
axs[1].plot(time_list,omy)
axs[1].set(ylabel='omega y')
axs[2].plot(time_list,omz)
axs[2].set(ylabel='omega z')
axs[2].set(xlabel='time')
#%%

nue_array=np.zeros([len(time_list)-1,10])
violation_list=[0]
for i in range(0,len(nue_list)):
    nue_array[i,:]=nue_list[i]
    violation_list.append(np.linalg.norm(nue_array[i,:]))

plt.figure()
plt.plot(time_list,violation_list)
plt.ylabel('2-Norm of Velocity Violation')
plt.xlabel('Time')
