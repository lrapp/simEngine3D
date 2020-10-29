import numpy as np
from constraint_value_functions import phi_dp1,phi_dp2,phi_cd,phi_d
from simEngine3D_dataload import data_file, DP1_PHI_partials,CD_PHI_partials,DP2_PHI_partials,D_PHI_partials, body
from nue_functions import DP1_nue,CD_nue
from gamma_functions import gamma_DP1,gamma_CD


#%%
def build_p_from_A(A):
    import numpy as np
    e0=np.sqrt((np.trace(A)+1)/4)
    if e0==0:
        print("e0=0")
    e1=(A[2,1]-A[1,2])/(4*e0)
    e2=(A[0,2]-A[2,0])/(4*e0)
    e3=(A[1,0]-A[0,1])/(4*e0)        

    p=np.empty([4,1])
    p[:,0]=np.array([e0,e1,e2,e3])

    return p
#%%


def build_p(theta):
    import numpy as np
    A=np.zeros([3,3])
    A[0,2]=-1
     
    A[1,0]=np.sin(theta)
    A[1,1]=-np.sin(theta)
    
    A[2,0]=-np.cos(theta)
    A[2,1]=-np.sin(theta)
    
    e0=np.sqrt((np.trace(A)+1)/4)
    if e0==0:
        print("e0=0")
    e1=(A[2,1]-A[1,2])/(4*e0)
    e2=(A[0,2]-A[2,0])/(4*e0)
    e3=(A[1,0]-A[0,1])/(4*e0)        

    p=np.empty([4,1])
    p[:,0]=np.array([e0,e1,e2,e3])

    return p


def build_A(p):
        import numpy as np
        A=np.empty([3,3])
           
        e0=p[0]
        e1=p[1]
        e2=p[2]
        e3=p[3]
        
        A[0,0]=2*((e0**2+e1**2)-0.5)
        A[0,1]=2*(e1*e2-e0*e3)
        A[0,2]=2*(e1*e3+e0*e2)
        
        A[1,0]=2*(e1*e2+e0*e3)
        A[1,1]=2*((e0**2+e2**2)-0.5)        
        A[1,2]=2*(e2*e3-e0*e1)
        
        A[2,0]=2*(e1*e3-e0*e2)
        A[2,1]=2*(e2*e3+e0*e1)
        A[2,2]=2*((e0**2+e3**2)-0.5)        
        
        return A
        

        
#%%    
def build_G(p):
        e0=p[0]
        e1=p[1]
        e2=p[2]
        e3=p[3]
        
        G=np.zeros([3,4])
        G[0,0]=-e1
        G[0,1]=e0
        G[0,2]=e3
        G[0,3]=-e2
        G[1,0]=-e2
        G[1,1]=-e3
        G[1,2]=e0
        G[1,3]=e1
        G[2,0]=-e3
        G[2,1]=e2
        G[2,2]=-e1
        G[2,3]=e0
        return G
    
def build_E(p):
        e0=p[0]
        e1=p[1]
        e2=p[2]
        e3=p[3]
        
        E=np.zeros([3,4])
        E[0,0]=-e1
        E[0,1]=e0
        E[0,2]=-e3
        E[0,3]=e2
        E[1,0]=-e2
        E[1,1]=e3
        E[1,2]=e0
        E[1,3]=-e1
        E[2,0]=-e3
        E[2,1]=-e2
        E[2,2]=e1
        E[2,3]=e0
        return E
    
 #%%
def build_B(p_vector,a_vector):
    e0=p_vector[0]
    e=p_vector[1:]
    a_bar=a_vector
    I3=np.identity(3)
    
    a_bar_T=np.transpose(a_bar)
    a_bar_tilde=tilde(a_bar)
    e_tilde=tilde(e)
    
    B=np.zeros([3,4])
    X=(2*(np.dot(np.dot(float(e0),I3)+e_tilde,a_bar)))
       
    Y=(2*(np.dot(e,a_bar_T)-np.dot(np.dot(float(e0),I3)+e_tilde,a_bar_tilde)))
    for i in range(0,len(X)):
        B[i,0]=X[i]
    B[0,1:]=Y[0]
    B[1,1:]=Y[1]
    B[2,1:]=Y[2]

    return B

 #%%
def tilde(a_vector):
    t=np.zeros([3,3],dtype=float)
    t[0,1]=-a_vector[2]
    t[0,2]=a_vector[1]
    t[1,0]=a_vector[2]
    t[1,2]=-a_vector[0]
    t[2,0]=-a_vector[1]
    t[2,1]=a_vector[0]
    return t
    #%%   

def build_A_angles(phi,theta,psi):
    import numpy as np
    A=np.zeros([3,3])    
    A[0,0]=np.cos(psi)*np.cos(phi)-np.cos(theta)*np.sin(phi)*np.sin(psi)
    A[0,1]=-np.sin(psi)*np.cos(phi)-np.cos(theta)*np.sin(phi)*np.cos(psi)
    A[0,2]=np.sin(theta)*np.sin(phi)
    
    A[1,0]=np.cos(psi)*np.sin(phi)+np.cos(theta)*np.cos(phi)*np.cos(psi)
    A[1,1]=-np.sin(psi)*np.sin(phi)+np.cos(theta)*np.cos(phi)*np.cos(psi)
    A[1,2]=-np.sin(theta)*np.cos(phi)
    
    A[2,0]=np.sin(theta)*np.sign(psi)
    A[2,1]=np.sin(theta)*np.cos(psi)
    A[2,2]=np.cos(theta)
    
    return A

def check_phi(X,cl,time,tol):
    phi=calc_phi(X,cl,time)
    check=0
    for i in phi:
        if i > tol:
            check=1
    if check == 1:
        print("initial conditions did not satisfy phi=0")
    else:
        print("initial conditions satisfy phi=0")
    return 
    


def calc_phi(X,cl,t):
    phi_values=[]
    for i in cl:
        if i.type.strip(" ' ") == 'CD':
            phi_values.append(phi_cd(X,i,eval(i.f)))
        elif i.type.strip(" ' ") == 'DP1':
            phi_values.append(phi_dp1(X,i,eval(i.f)))
    
        #Add euler parameter constraint
    for k in range(1,len(X)):
        phi_values.append(round(float((np.dot(np.transpose(X[k].q[3:]),X[k].q[3:])-1)),12))
    return phi_values

def calc_partials(X,cl):
    phi_partials_values=[]
    for i in cl:
            if i.type.strip(" ' ") == 'CD':
                dri,dpi,drj,dpj=CD_PHI_partials(X,i)
                phi_partials_values.append([drj,dpj])
            if i.type.strip(" ' ") == 'DP1':
                dri,dpi,drj,dpj=DP1_PHI_partials(X,i)
                phi_partials_values.append([drj,dpj])    
    return phi_partials_values


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
def calc_phi_HW82(X,cl,t):
    phi_values=[]
    for i in cl:
        if i.type.strip(" ' ") == 'CD':
            phi_values.append(phi_cd(X,i,eval(i.f)))
        elif i.type.strip(" ' ") == 'DP1':
            phi_values.append(phi_dp1(X,i,eval(i.f)))
    return phi_values
#%%

def calc_nue(X,cl,df):
    nue_values=[]
    for i in cl:
            if i.type.strip(" ' ") == 'CD':
                nue_values.append(CD_nue(X,i,df))
            if i.type.strip(" ' ") == 'DP1':
                nue_values.append(DP1_nue(X,i,df))            
    return nue_values

#%%

def calc_gamma(X,cl,ddf):
    gamma_values=[]
    for i in cl:
            if i.type.strip(" ' ") == 'CD':
                gamma_values.append(gamma_CD(X,i,ddf))
            if i.type.strip(" ' ") == 'DP1':
                gamma_values.append(gamma_DP1(X,i,ddf))
    return gamma_values

#%%
def build_ja(X,phi_partials_values):
    ja_list=[]
    for i in range(0,len(phi_partials_values)):
        ca=np.zeros([7])
        #order for [d_ri d_pi d_rj d_pj]
        ca[0:3]=phi_partials_values[i][0]
        ca[3:7]=phi_partials_values[i][1]
    
        ca.shape=(1,7)
        ja_list.append(ca)
    
    
    jacobian=np.zeros([len(phi_partials_values),7])
    for i in range(0,len(ja_list)):
        jacobian[i,:]=ja_list[i]
        
    #add euler parameter normalization constraint to jacobian
#    for k in range(1,len(X)):
#        jacobian[5+k,3:]=2*X[k].q[3:].T
                
    return jacobian