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

# Read in body and constraint information
file="C:\\Users\\Logan\\Desktop\\simEngine3D\\HW8_P2_input.txt"
X=data_file(file)
constraint_list=constraints_in(file)

#%%
h=0.001
L=2
xa=0.05*0.05
volume=xa*L
rho=7800
m1=rho*volume
m2=rho*volume/2

counter=0
time_list=np.arange(0,5+h,h)
                       
time=time_list[0]

#df=np.pi*np.sin(2*time)*np.cos((np.pi*np.cos(time*2))/4)*-1/2
#ddf=np.pi**2*np.sin(time*2)**2*np.sin((np.pi*np.cos(time*2))/4)*(-1/4)-np.pi*np.cos(time*2)*np.cos((np.pi*np.cos(time*2))/4)
df=0
ddf=0

nb=3
nc=len(constraint_list)

r_list=np.zeros([3*nb,len(time_list)])
r_d_list=np.zeros([3*nb,len(time_list)])
r_dd_list=[]

p_list=np.zeros([4*nb,len(time_list)])
p_d_list=np.zeros([4*nb,len(time_list)])
p_dd_list=[]


lambda_p_list=[]
lagrange_list=[]
         
F0=np.array([0,0,0]) 
F0.shape=(3,1)
F1=np.array([0,0,m1*-9.81]) 
F1.shape=(3,1)
F2=np.array([0,0,m2*-9.81]) 
F2.shape=(3,1)

#F=np.zeros([3*nb,1])
#F[0:3]=F1
#F[3:6]=F2
 
F=np.zeros([3*nb,1])
F[0:3]=F0
F[3:6]=F1
F[6:9]=F2

tau=np.array([0,0,0])
tau.shape=(3,1)

#%% set body j starting values
r_start=np.array([0,L,0]) 
A_start=np.zeros([3,3])
A_start[0,2]=np.cos(0)
A_start[1,0]=np.cos(0)
A_start[2,1]=np.cos(0)
p_start=build_p_from_A(A_start)


x0_q_start=np.array([0,0,0,1,0,0,0])
x0_q_start.shape=(7,1)
X[0].q=x0_q_start


r_start.shape=(3,1)
p_start.shape=(4,1)
X[1].q[0:3]=r_start
X[1].q[3:]=p_start 
X[1].A_rotation=build_A(X[1].q[3:])

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

#%%
tol=1e-5
phi=calc_phi_HW82(X,constraint_list,0)
check_phi(X,constraint_list,0,tol)
PHI_P=np.zeros([nb,1])
for i in range(0,len(X)):
    PHI_P[i-1]=1/2*(X[i].q[3:].T @ X[i].q[3:]) - 1/2
         
if PHI_P.any() > tol:
    print("initial conditions do not satisfy PHI_P=0")
else:
    print("initial conditions satisfy PHI_P=0")


r_dd_list=[]
p_dd_list=[]
lambda_p_list=[]
lagrange_list=[]
#%%

partials=calc_partials_HW82(X,constraint_list)
jacobian=np.vstack(partials)


nue_values=calc_nue(X,constraint_list,df)


phi_r=np.zeros([nc,3*nb])
#phi_r=jacobian[:,3:9]
phi_r=jacobian[:,0:3*nb]

phi_p=np.zeros([nc,4*nb])
phi_p=jacobian[:,3*nb:3*nb+4*nb]
#phi_p=jacobian[:,13:]

gamma_values=calc_gamma(X,constraint_list,ddf)
           
M0=np.zeros([3,3])
M1=m1*np.identity(3)
M2=m2*np.identity(3)
M_list=[M0,M1,M2]

M=np.zeros([3*nb,3*nb])
for k in range(0,nb):
    M[k*nb:k*nb+3,k*nb:k*nb+3]=M_list[k]

 
J_bar1=np.zeros([3,3])

b=0.05/2
c=0.05/2

J_bar0=np.zeros([3,3])

J_bar1[0,0]=1/12*m1*(b**2+c**2)
J_bar1[1,1]=1/12*m1*(L**2+c**2)
J_bar1[2,2]=1/12*m1*(L**2+b**2)

J_bar2=np.zeros([3,3])
J_bar2[0,0]=1/12*m1*(b**2+c**2)
J_bar2[1,1]=1/12*m1*(L/2**2+c**2)
J_bar2[2,2]=1/12*m1*(L/2**2+b**2)

J_bar=[J_bar0,J_bar1,J_bar2]

G=[]
for i in range(0,len(X)):
    G.append(build_G(X[i].q[3:]))

J_P_list=[]
for i in range(0,len(G)):
    J_P_list.append(4*np.dot(np.dot(np.transpose(G[i]),J_bar[i]),G[i]))    
    
J_P=np.zeros([nb*4,nb*4])
for k in range(0,nb):
    J_P[k*(4):k*(4)+4,k*(4):k*(4)+4]=J_P_list[k]
    

G_dot=[]
for i in range(0,len(X)):
    G_dot.append(build_G(X[i].p_dot))

tau_hat_list=[]
for i in range(0,len(G_dot)):
    tau_hat_list.append(8*np.dot(np.dot(np.dot(G_dot[i].T,J_bar[i]),G_dot[i]),X[i].q[3:]))

tau_hat=np.zeros([4*nb,1])
for k in range(0,nb):
    tau_hat[k*4:k*4+4]=tau_hat_list[k]

P=np.zeros([nb,nb*4])
for k in range(0,nb):
    P[k,k*4:k*4+4]=X[k].q[3:].T


LHS=np.zeros([8*nb+nc,8*nb+nc])
LHS[0:3*nb,0:3*nb]=M
LHS[0:3*nb,8*nb:8*nb+nc]=phi_r.T 
   
LHS[3*nb:7*nb,3*nb:7*nb]=J_P
LHS[3*nb:7*nb,7*nb:7*nb+nb]=P.T
LHS[3*nb:7*nb,8*nb:8*nb+nc]=phi_p.T

LHS[7*nb:8*nb,3*nb:7*nb]=P
   
LHS[8*nb:8*nb+nc,0:3*nb]=phi_r
LHS[8*nb:8*nb+nc,3*nb:7*nb]=phi_p
   
RHS=np.zeros([8*nb+nc,1])

gamma_p=np.zeros([nb,1])
for i in range(0,len(X)):   
    gamma_p[i]=(-2*np.dot(X[i].p_dot.T,X[i].p_dot)) #slide 473

RHS[0:3*nb]=F
RHS[3*nb:7*nb]=tau_hat
RHS[7*nb:8*nb]=gamma_p
gamma_hat=np.array(gamma_values)
gamma_hat.shape=(10,1)
RHS[8*nb:8*nb+nc]=gamma_hat
   
unknowns=np.dot(np.linalg.inv(LHS),RHS) 
r_dd_list.append(unknowns[0:6])
p_dd_list.append(unknowns[6:14])
lambda_p_list.append(unknowns[14:16])
lagrange_list.append(unknowns[16:26])
  


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


#%%
n=1
time=time_list[n]

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
    
    #%%to do: loop this
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
    
    nue_values=calc_nue(X_n,constraint_list,df)
    
    phi_r=np.zeros([10,6])
    phi_r[:,0:3]=jacobian[:,0:3] #get ri
    phi_r[:,3:6]=jacobian[:,3:6]
    
    phi_p=np.zeros([10,8])
    phi_p[:,0:4]=jacobian[:,6:10]
    phi_p[:,4:8]=jacobian[:,10:14]    
    
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
    P[0,0:4]=X_n[1].q[3:].T
    P[1,4:8]=X_n[2].q[3:].T
            
     
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
    
    delta_z=np.dot(np.linalg.inv(PSI),-g)           
    z=z+delta_z
       
       
