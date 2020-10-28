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

#%% Read in body and constraint information
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


r_list=np.zeros([3,len(time_list)])
r_d_list=np.zeros([3,len(time_list)])
r_dd_list=[]

p_list=np.zeros([4,len(time_list)])
p_d_list=np.zeros([4,len(time_list)])
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

#%% set body j starting values
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
PHI_P=np.zeros([2,1])
for i in range(1,len(X)):
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
#jacobian=build_ja(X,partials)

nue_values=calc_nue(X,constraint_list,df)

#jacobian=jacobian[:-2,:] #remove euler parameterization constraint

phi_r=np.zeros([10,6])
phi_r[:,0:3]=jacobian[:,0:3] #get ri
phi_r[:,3:6]=jacobian[:,3:6]

phi_p=np.zeros([10,8])
phi_p[:,0:4]=jacobian[:,6:10]
phi_p[:,4:8]=jacobian[:,10:14]

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
r_dd_list.append(unknowns[0:3])
p_dd_list.append(unknowns[6:10])
lambda_p_list.append(unknowns[14:16])
lagrange_list.append(unknowns[16:26])
  
