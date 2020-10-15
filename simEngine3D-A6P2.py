#simEngine3D-A6P2
import os
os.chdir("C:\\Users\\Logan\\Desktop\\simEngine3D")
from constraint_value_functions import phi_dp1,phi_dp2,phi_cd,phi_d
from nue_functions import DP1_nue,CD_nue,DP2_nue,D_nue
from gamma_functions import gamma_DP1,gamma_CD,gamma_DP2, gamma_D
from simEngine3D_dataload import data_file, DP1_PHI_partials,CD_PHI_partials,DP2_PHI_partials,D_PHI_partials, body
from HW6_Q2_flags import flags_in,points_in, f_in

from constraints_in import constraints_in
from simEngine3D_functions import build_p


constraint_list=constraints_in()

X=data_file()
import numpy as np

#file="C:\\Users\\Logan\\Desktop\\ME-751\\HW6\\input.txt" 

theta=0
L=2
t=0
#%%
theta=np.pi/4

p_j=build_p(theta)

#%%
for i in range(0,len(X)):
    if X[i].name == "body j":
        X[i].q[3:]=p_j
           
         
#X[1].q[:3,0]=np.array([0,L*np.sin(theta),-L*np.cos(theta)])

X[1].q[:3,0]=np.array([0,1.2,-1.2])


ft=np.pi/4*np.cos(2*t)
df=-np.pi/2*np.sin(2*t)
ddf=-np.pi*np.cos(2*t)

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
            phi_partials_values.append(CD_PHI_partials(X,i))
        if i.type.strip(" ' ") == 'DP1':
            phi_partials_values.append(DP1_PHI_partials(X,i))
            

           
nue_values=[]
for i in constraint_list:
        if i.type.strip(" ' ") == 'CD':
            nue_values.append(CD_nue(X,i,0))
        if i.type.strip(" ' ") == 'DP1':
            nue_values.append(DP1_nue(X,i,0))            
            
gamma_values=[]
for i in constraint_list:
        if i.type.strip(" ' ") == 'CD':
            gamma_values.append(gamma_CD(X,i,0))
        if i.type.strip(" ' ") == 'DP1':
            gamma_values.append(gamma_DP1(X,i,0))


print("PHI(q,t)=",str(phi_values))

for i in range(0,len(phi_partials_values)):
    print("phi_qi_cd_",str(constraint_list[i].ID),"=",str(phi_partials_values[i][0]),str(","),str(phi_partials_values[i][1]))
    print("phi_qj_cd_",str(constraint_list[i].ID),"=",str(phi_partials_values[i][2]),str(","),str(phi_partials_values[i][3]))
    
    
#%%
ja_list=[]
for i in phi_partials_values:
    ca=np.zeros([14])
    ca[0:3]=i[0]
    ca[3:7]=i[1]
    ca[7:10]=i[2]
    ca[10:14]=i[3]
    ja_list.append(ca)


jacobian=np.zeros([7,14])
for i in range(0,len(ja_list)):
    jacobian[i,:]=ja_list[i]
    
#add euler parameter normalization constraint to jacobian
p_i=X[0].q[3:]
jacobian[6,3:7]=2*np.transpose(p_i)
jacobian[6,10:]=2*np.transpose(p_j)

#%%
#Newton Raphson method
q0=X[1].q
phi=np.array(phi_values)    

q_new=q0-np.linalg.lstsq(jacobian,phi)






#%%
   
print("nue=",str(nue_values))
print("gamma=",str(gamma_values))

