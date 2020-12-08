#%%
# calculate p given A
#def get_p(A):
#    e0 = np.sqrt(0.25 * (np.trace(A) + 1.))
#    e1 = np.sqrt(0.25 * (-np.trace(A) + 1. + 2.* A[0][0]))
#    e2 = np.sqrt(0.25 * (-np.trace(A) + 1. + 2.* A[1][1]))
#    e3 = np.sqrt(0.25 * (-np.trace(A) + 1. + 2.* A[2][2]))
#    return np.array([e0, e1, e2, e3])

#%%
#def build_p(theta):
#    import numpy as np
#    A=np.zeros([3,3])
#    A[0,2]=-1
#     
#    A[1,0]=np.sin(theta)
#    A[1,1]=-np.sin(theta)
#    
#    A[2,0]=-np.cos(theta)
#    A[2,1]=-np.sin(theta)
#    
#    e0=np.sqrt((np.trace(A)+1)/4)
#    if e0==0:
#        print("e0=0")
#    e1=(A[2,1]-A[1,2])/(4*e0)
#    e2=(A[0,2]-A[2,0])/(4*e0)
#    e3=(A[1,0]-A[0,1])/(4*e0)        
#
#    p=np.empty([4,1])
#    p[:,0]=np.array([e0,e1,e2,e3])
#
#    return p
#%%
#
#def build_A_angles(phi,theta,psi):
#    import numpy as np
#    A=np.zeros([3,3])    
#    A[0,0]=np.cos(psi)*np.cos(phi)-np.cos(theta)*np.sin(phi)*np.sin(psi)
#    A[0,1]=-np.sin(psi)*np.cos(phi)-np.cos(theta)*np.sin(phi)*np.cos(psi)
#    A[0,2]=np.sin(theta)*np.sin(phi)
#    
#    A[1,0]=np.cos(psi)*np.sin(phi)+np.cos(theta)*np.cos(phi)*np.cos(psi)
#    A[1,1]=-np.sin(psi)*np.sin(phi)+np.cos(theta)*np.cos(phi)*np.cos(psi)
#    A[1,2]=-np.sin(theta)*np.cos(phi)
#    
#    A[2,0]=np.sin(theta)*np.sign(psi)
#    A[2,1]=np.sin(theta)*np.cos(psi)
#    A[2,2]=np.cos(theta)
#    
#    return A
#%%
#
#def calc_partials2(X,cl):
#    if X[1].ground:
#        r_cols=[0,3]
#        p_cols=[3,7]
#        cols=7
#    else:
#        cols=7*len(X)
#        
#    phi_partials_values=np.zeros([len(cl),cols])
#    counter=0
#    for i in cl:
#            if i.type.strip(" ' ") == 'CD':
#                phi_partials_values[counter,r_cols[0]:r_cols[1]]=CD_PHI_r(X,i)
#                phi_partials_values[counter,p_cols[0]:p_cols[1]]=CD_PHI_p(X,i)
#                
#            if i.type.strip(" ' ") == 'DP1':
#
#                
#    return phi_partials_values
#%%

def calc_partials_HW82(X,cl):
    phi_partials_values=np.zeros([len(cl),7*len(X)])
    counter=0
    for i in cl:
            if i.type.strip(" ' ") == 'CD':
                dri,dpi,drj,dpj=CD_PHI_partials(X,i)
                if i.body_i_name == "body i":
                    phi_partials_values[counter,0:3]=dri
                    phi_partials_values[counter,3*len(X):3*len(X)+4]=dpi                          
                                       
                if i.body_i_name == "body j":
                    phi_partials_values[counter,3:6]=dri          
                    phi_partials_values[counter,3*len(X)+4:3*len(X)+8]=dpi
                                        
                if i.body_j_name == "body j":
                    phi_partials_values[counter,3:6]=drj          
                    phi_partials_values[counter,3*len(X)+4:3*len(X)+8]=dpj       
                                        
                if i.body_j_name == "body k":
                    phi_partials_values[counter,6:9]=drj          
                    phi_partials_values[counter,3*len(X)+8:3*len(X)+12]=dpj                  

            if i.type.strip(" ' ") == 'DP1':
                
                dri,dpi,drj,dpj=DP1_PHI_partials(X,i)

                if i.body_i_name == "body i":
                    phi_partials_values[counter,0:3]=dri
                    phi_partials_values[counter,3*len(X):3*len(X)+4]=dpi                          
                                       
                if i.body_i_name == "body j":
                    phi_partials_values[counter,3:6]=dri          
                    phi_partials_values[counter,3*len(X)+4:3*len(X)+8]=dpi
                                        
                if i.body_j_name == "body j":
                    phi_partials_values[counter,3:6]=drj          
                    phi_partials_values[counter,3*len(X)+4:3*len(X)+8]=dpj       
                                        
                if i.body_j_name == "body k":
                    phi_partials_values[counter,6:9]=drj          
                    phi_partials_values[counter,3*len(X)+8:3*len(X)+12]=dpj                                                    
        
            counter=counter+1
                
    return phi_partials_values

#%%
#%%
def calc_phi_HW82(X,cl,t):
    phi_values=[]
    for i in cl:
        if i.type.strip(" ' ") == 'CD':
            phi_values.append(phi_cd(X,i,eval(i.f)))
        elif i.type.strip(" ' ") == 'DP1':
            phi_values.append(phi_dp1(X,i,eval(i.f)))
    return phi_values