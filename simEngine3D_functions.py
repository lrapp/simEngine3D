import numpy as np
from constraint_value_functions import phi_dp1,phi_dp2,phi_cd,phi_d
from simEngine3D_dataload import data_file, DP1_PHI_partials,CD_PHI_partials,DP2_PHI_partials,D_PHI_partials, body


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


def calc_phi(X,cl,t):
    phi_values=[]
    for i in cl:
        if i.type.strip(" ' ") == 'CD':
            phi_values.append(phi_cd(X,i,eval(i.f)))
        elif i.type.strip(" ' ") == 'DP1':
            phi_values.append(phi_dp1(X,i,eval(i.f)))
    
        #Add euler parameter constraint
    phi_values.append(round(float((np.dot(np.transpose(X[1].q[3:]),X[1].q[3:])-1)),12))
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

def build_ja(X,phi_partials_values):
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
    return jacobian