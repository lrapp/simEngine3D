#simEngine3D-A6P2
import os
os.chdir("C:\\Users\\Logan\\Desktop\\simEngine3D")
from constraint_value_functions import phi_dp1,phi_dp2,phi_cd,phi_d
from nue_functions import DP1_nue,CD_nue,DP2_nue,D_nue
from gamma_functions import gamma_DP1,gamma_CD,gamma_DP2, gamma_D
from simEngine3D_dataload import data_file, DP1_PHI_partials,CD_PHI_partials,DP2_PHI_partials,D_PHI_partials, body
from HW6_Q2_flags import flags_in,points_in, f_in

from constraints_in import constraints_in
from simEngine3D_functions import build_p, build_A, calc_phi, calc_partials, build_ja


constraint_list=constraints_in()

X=data_file()
import numpy as np

#file="C:\\Users\\Logan\\Desktop\\ME-751\\HW6\\input.txt" 


L=2
t=0
#%%
theta=np.pi/4

#


p_j=build_p(theta)
X[1].A_rotation=build_A(p_j)

df=np.cos(-1/2**np.pi*np.sin(2*t))
ddf=np.cos(-np.pi*np.cos(2*t))
#%%
for i in range(0,len(X)):
    if X[i].name == "body j":
        X[i].q[3:]=p_j
           
         
X[1].q[:3,0]=np.array([0,L*np.sin(theta),-L*np.cos(theta)])


phi_values=[]
for i in constraint_list:
    if i.type.strip(" ' ") == 'CD':
        phi_values.append(phi_cd(X,i,eval(i.f)))
    elif i.type.strip(" ' ") == 'DP1':
        phi_values.append(phi_dp1(X,i,eval(i.f)))

#Add euler parameter constraint
phi_values.append(round(float((np.dot(np.transpose(X[1].q[3:]),X[1].q[3:])-1)),12))

phi_partials_values=[]
for i in constraint_list:
        if i.type.strip(" ' ") == 'CD':
            dri,dpi,drj,dpj=CD_PHI_partials(X,i)
            phi_partials_values.append([drj,dpj])
        if i.type.strip(" ' ") == 'DP1':
            dri,dpi,drj,dpj=DP1_PHI_partials(X,i)
            phi_partials_values.append([drj,dpj])
            

           
nue_values=[]
for i in constraint_list:
        if i.type.strip(" ' ") == 'CD':
            nue_values.append(CD_nue(X,i,df))
        if i.type.strip(" ' ") == 'DP1':
            nue_values.append(DP1_nue(X,i,df))            
            
gamma_values=[]
for i in constraint_list:
        if i.type.strip(" ' ") == 'CD':
            gamma_values.append(gamma_CD(X,i,ddf))
        if i.type.strip(" ' ") == 'DP1':
            gamma_values.append(gamma_DP1(X,i,ddf))


print("PHI(q,t)=",str(phi_values))

#for i in range(0,len(phi_partials_values)):
#    print("phi_qi_cd_",str(constraint_list[i].ID),"=",str(phi_partials_values[i][0]),str(","),str(phi_partials_values[i][1]))
#    print("phi_qj_cd_",str(constraint_list[i].ID),"=",str(phi_partials_values[i][2]),str(","),str(phi_partials_values[i][3]))
    
    
#%% Assembel Jacobian
ja_list=[]
for i in range(0,len(phi_partials_values)):
    ca=np.zeros([7])
    #order for [d_ri d_pi d_rj d_pj]
    ca[0:3]=phi_partials_values[i][0]
    ca[3:7]=phi_partials_values[i][1]

    ca.shape=(1,7)
    ja_list.append(ca)


jacobian=np.zeros([7,7])
for i in range(0,len(ja_list)):
    jacobian[i,:]=ja_list[i]
    
#add euler parameter normalization constraint to jacobian
jacobian[6,3:]=2*X[1].q[3:].T
#%%
#Newton Raphson method
q0=np.zeros([7,1])
q0[:]=X[1].q

phi=np.array(phi_values)    
phi.shape=(7,1)

error=10
tol=1e-5
counter=1

while abs(error) > tol:
    q_new=q0-np.linalg.lstsq(jacobian,phi)[0]
    error = np.linalg.norm(q0-q_new)
    
    phi_new=[]
    X[1].q=q_new #Assing new q to body J
    X[1].A_rotation=build_A(q_new[3:])     #Build new rotation matrix from new p values
    phi=np.array(calc_phi(X,constraint_list,t)) #calculate new values of phi
    phi.shape=(7,1)
    partials_new=calc_partials(X,constraint_list) #calculate new partials
    jacobian=build_ja(X,partials_new) #build new jacobian
    q0=q_new
    counter=counter+1
    



#%%
print("nue=",str(nue_values))
print("gamma=",str(gamma_values))


