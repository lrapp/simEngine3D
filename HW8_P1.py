import os
os.chdir("C:\\Users\\Logan\\Desktop\\simEngine3D")
from pendulum_function import pendulum
import numpy as np
from constraints_in import constraints_in
from simEngine3D_dataload import data_file, DP1_PHI_partials,CD_PHI_partials,DP2_PHI_partials,D_PHI_partials, body
from simEngine3D_functions import build_p, build_A, calc_phi, calc_partials, build_ja,build_G,tilde,calc_nue,check_phi,calc_gamma,build_E
import matplotlib.pyplot as plt
from new_partials import DP1_phi_parital_lagrange,CD_phi_parital_lagrange
import time as ttime

tic=ttime.perf_counter()
h=0.001


L=2
xa=0.05*0.05
volume=xa*L
rho=7800
m=rho*volume

counter=0
time_list=np.arange(0,5+h,h)
                       
time=time_list[0]

df=np.pi*np.sin(2*time)*np.cos((np.pi*np.cos(time*2))/4)*-1/2
ddf=np.pi**2*np.sin(time*2)**2*np.sin((np.pi*np.cos(time*2))/4)*(-1/4)-np.pi*np.cos(time*2)*np.cos((np.pi*np.cos(time*2))/4)

r_list=np.zeros([3,len(time_list)])
r_d_list=np.zeros([3,len(time_list)])

r_dd_list=[]


p_list=np.zeros([4,len(time_list)])
#p_d_list=[]
p_d_list=np.zeros([4,len(time_list)])

p_dd_list=[]

lambda_p_list=[]
lagrange_list=[]
         

X=data_file()
constraint_list=constraints_in()

F=np.array([0,0,m*-9.81]) 
F.shape=(3,1)

tau=np.array([0,0,0])
tau.shape=(3,1)


r_start=np.array([0,1.41421356,-1.41421356]) 
p_start=np.array([.653281459,.270598060,.653281506,.27059804])
r_start.shape=(3,1)
p_start.shape=(4,1)
X[1].q[0:3]=r_start #set initial position from HW7 
X[1].q[3:]=p_start #set initial position from HW7 
X[1].A_rotation=build_A(X[1].q[3:])


#%%
#check initial conditions
tol=1e-5
check_phi(X,constraint_list,0,tol)
PHI_P=1/2*(X[1].q[3:].T @ X[1].q[3:]) - 1/2
if PHI_P > tol:
    print("initial conditions do not satisfy PHI_P=0")
else:
    print("initial conditions satisfy PHI_P=0")

partials=calc_partials(X,constraint_list)    
jacobian=build_ja(X,partials)
nue_values=calc_nue(X,constraint_list,df)

phi_r=jacobian[:,0:3]
phi_p=jacobian[:,3:]

phi_r=phi_r[0:6]
phi_p=phi_p[0:6]

gamma_values=calc_gamma(X,constraint_list,ddf)

M=m*np.identity(3)
J_bar=np.zeros([3,3])

b=0.05/2
c=0.05/2

J_bar[0,0]=1/12*m*(b**2+c**2)
J_bar[1,1]=1/12*m*(L**2+c**2)
J_bar[2,2]=1/12*m*(L**2+b**2)

G=build_G(X[1].q[3:])
           
J_P=4*np.dot(np.dot(np.transpose(G),J_bar),G)    
P=np.zeros([1,4])
P[0,:]=np.transpose(X[1].q[3:])

#construct blue LHS matrix on slide 640 of 2019 slides
LHS=np.zeros([14,14])
LHS[0:3,0:3]=M
LHS[0:3,8:15]=phi_r.T
LHS[3:7,3:7]=J_P
LHS[3:7,7:8]=P.T
LHS[3:7,8:15]=phi_p.T
LHS[7:8,3:7]=P
LHS[8:15,0:3]=phi_r
LHS[8:15,3:7]=phi_p
   
   
G_dot=build_G(X[1].p_dot)
tau_hat=8*np.dot(np.dot(np.dot(G_dot.T,J_bar),G_dot),X[1].q[3:])
gamma_p=-2*np.dot(X[1].p_dot.T,X[1].p_dot) #slide 473

gamma_hat=np.array(gamma_values[0:6])
gamma_hat.shape=(6,1)
                 
RHS=np.zeros([14,1])
RHS[0:3]=F
RHS[3:7]=tau_hat
RHS[7:8]=gamma_p
RHS[8:14]=gamma_hat


unknowns=np.dot(np.linalg.inv(LHS),RHS)    
r_dd_list.append(unknowns[0:3])
p_dd_list.append(unknowns[3:7])
lambda_p_list.append(unknowns[7:8])
lagrange_list.append(unknowns[8:14])


r_list[:,0:1]=X[1].q[0:3]

p_list[:,0:1]=X[1].q[3:]

r_d_start=X[1].r_dot
r_d_start.shape=(3,1)
r_d_list[:,0:1]=r_d_start
        
#p_d_list.append(X[1].p_dot)
p_d_start=X[1].p_dot
p_d_start.shape=(4,1)
p_d_list[:,0:1]=p_d_start

#%%
#start steps on slide 641
n=1
time=time_list[n]
for ii in range(1,len(time_list)):
    
    time=time_list[ii]
    df=np.pi*np.sin(2*time)*np.cos((np.pi*np.cos(time*2))/4)*-1/2
    ddf=np.pi**2*np.sin(time*2)**2*np.sin((np.pi*np.cos(time*2))/4)*(-1/4)-np.pi*np.cos(time*2)*np.cos((np.pi*np.cos(time*2))/4)

    r_dd=r_dd_list[n-1]
    p_dd=p_dd_list[n-1]
    lagrange=lagrange_list[n-1]
    lambda_p=lambda_p_list[n-1]
    
    z_0=np.zeros([14,1])
    z_0[0:3]=r_dd_list[n-1]
    z_0[3:7]=p_dd_list[n-1]
    z_0[7:8]=lambda_p_list[n-1]
    z_0[8:14]=lagrange_list[n-1]
    
    tol=1e-4
    error=1
    count=1
    while abs(error) > tol and count < 30:
        
        if count == 1:
            if n==1: #use first order to seed 
                z=z_0
                beta_0=1
                alpha1=1
                C_r_n=alpha1*r_list[:,n-1]
                C_r_n.shape=(3,1)
                C_p_n=alpha1*p_list[:,n-1]
                C_p_n.shape=(4,1)
                
                C_r_dot_n=alpha1*r_d_list[:,n-1]
                C_r_dot_n.shape=(3,1)
                C_p_dot_n=alpha1*p_d_list[:,n-1]
                C_p_dot_n.shape=(4,1)                
            else:
                z=z_0
                beta_0=2/3
                C_r_n=4/3*r_list[:,n-1]-1/3*r_list[:,n-2]
                C_r_n.shape=(3,1)                
                C_p_n=4/3*p_list[:,n-1]-1/3*p_list[:,n-2]
                C_p_n.shape=(4,1)
                
                C_r_dot_n=4/3*r_d_list[:,n-1]-1/3*r_d_list[:,n-2]
                C_r_dot_n.shape=(3,1)
                C_p_dot_n=4/3*p_d_list[:,n-1]-1/3*p_d_list[:,n-2]
                C_p_dot_n.shape=(4,1)
    
                
        else:
            z=z
         
            
            
        lambda_p=z[7:8]             
        lagrange=z[8:14]
        
        r=C_r_n+(beta_0**2)*(h**2)*z[0:3]
        p=C_p_n+(beta_0**2)*(h**2)*z[3:7]        
    
        r_d=C_r_dot_n+beta_0*h*z[0:3]
        p_d=C_p_dot_n+beta_0*h*z[3:7]        
    
    
        #update body object with new values 
        X_n=X
        X_n[1].q[0:3]=r
        X_n[1].q[3:]=p
        X_n[1].A_rotation=build_A(X_n[1].q[3:])
        X_n[1].r_dot=r_d
        X_n[1].p_dot=p_d
        X_n[1].r_d_dot=z[0:3]
        X_n[1].p_d_dot=z[3:7]
         
        PHI=calc_phi(X_n,constraint_list,time)
        PHI=np.array(PHI[0:6]) #remove euler parameterization constraint
        
        partials=calc_partials(X_n,constraint_list)            
        jacobian_n=build_ja(X_n,partials)
        phi_r=jacobian_n[:,0:3]
        phi_p=jacobian_n[:,3:]
        
        phi_r=phi_r[0:6] #don't consider euler parameter contstraints
        phi_p=phi_p[0:6] #don't consider euler parameter constraints
                    
        PHI_P=1/2*(p.T @ p) - 1/2
        G_dot=build_G(p_d)
        tau_hat=8*np.dot(np.dot(np.dot(G_dot.T,J_bar),G_dot),p)
        
        
        G=build_G(X_n[1].q[3:])
        J_P=4*np.dot(np.dot(np.transpose(G),J_bar),G)    
        P=np.zeros([1,4])
        P[0,:]=np.transpose(X_n[1].q[3:])
        
        g=np.zeros([14,1])
        g[0:3]=(M @ z[0:3])+(phi_r.T @ lagrange)-F
        g[3:7]=J_P @ z[3:7] + phi_p.T @ lagrange + P.T @ lambda_p - tau_hat
        g[7:8]=1/((beta_0**2)*(h**2))*PHI_P
        last_g=(1/((beta_0**2)*(h**2)) * PHI)
        last_g.shape=(6,1)
        g[8:14]=last_g
    
          
         
        PSI=np.zeros([14,14])
        PSI[0:3,0:3]=M
        PSI[0:3,8:14]=phi_r.T
        
        PSI[3:7,4:8]=J_P
        PSI[3:7,7:8]=P.T
        PSI[3:7,8:14]=phi_p.T
        PSI[7:8,3:7]=P
        PSI[8:14,0:3]=phi_r
        PSI[8:14,3:7]=phi_p
        
    
        delta_z=-np.dot(np.linalg.inv(PSI),g)    
        z=z+delta_z
        error=np.linalg.norm(delta_z)
        count=count+1
        if count > 29:
            print("did not converge")
    
    r_list[:,n:n+1]=r
    p_list[:,n:n+1]=p

    r_d_list[:,n:n+1]=r_d
    p_d_list[:,n:n+1]=p_d
            
    r_dd_list.append(z[0:3])
    p_dd_list.append(z[3:7])
    lagrange_list.append(z[8:14])
    lambda_p_list.append(z[7:8])
    n=n+1

toc=ttime.perf_counter()
    
    #%%

x=r_list[0,:]
y=r_list[1,:]
z=r_list[2,:]

plt.figure()
plt.plot(time_list,x)


plt.figure()
plt.plot(time_list,y)

plt.figure()
plt.plot(time_list,z)

#%%calculate omega
E_list=[]
for i in range(0,max(p_list.shape)):
    E_list.append(build_E(p_list[:,i]))

omega_list=np.zeros([3,max(p_list.shape)])
for i in range(0,max(p_d_list.shape)):
    omega=2*(np.dot(E_list[i], p_d_list[:,i]))
    omega.shape=(3,1)
    omega_list[:,i:i+1]=omega

#%%
plt.plot(time_list,omega_list[0,:])
plt.plot(time_list,omega_list[1,:])
plt.plot(time_list,omega_list[2,:])

elapsed_time=toc-tic
print(elapsed_time)
