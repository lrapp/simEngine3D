import numpy as np
from simEngine3D_functions import tilde

#%%
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
def DP2_phi(X,C,ft):
    body_1=C.body_i_name
    body_2=C.body_j_name
    
    for i in X:
        if i.name == body_1:
            i_index=X.index(i)
        if i.name == body_2:
            j_index=X.index(i) 
               
    I=X[i_index]
    J=X[j_index]    
    
    a_bar_i=I.a_bar

    A_i=build_A(I.q[3:])
    A_j=build_A(J.q[3:])

    r_j=J.q[0:3]
    s_bar_j=J.s_bar
    
    r_i=I.q[0:3]
    s_bar_i=I.s_bar
    
    a_bar_i_T=np.transpose(a_bar_i)
    
    phi_dp2=np.dot(np.dot(a_bar_i_T,A_i),(r_j+np.dot(A_j,s_bar_j)-r_i-np.dot(A_i,s_bar_i)))
             

    return float(phi_dp2)

#%%

def DP2_nue(X,C,dft):

    if C.name == "driving":
        dft=dft
    else:
        dft=0
        
    return dft
#%%    
def gamma_DP2(X,ddf):
        #Set bodys i and j 
    for i in X:
        if i.name == "body i":
            i_index=X.index(i)
        if i.name == "body j":
            j_index=X.index(i) 
            
    I=X[i_index]
    J=X[j_index]  
    
    #Set a and A for bodies i and j
    a_bar_i=I.a_bar
    a_bar_j=J.a_bar
    A_i=build_A(I.q[3:])
    A_j=build_A(J.q[3:])
    
    p_i=I.q[3:]
    p_j=J.q[3:]
    p_dot_i=I.p_dot
    p_dot_j=J.p_dot

    a_i=np.dot(A_i,a_bar_i)
    a_j=np.dot(A_j,a_bar_j)
    
    
    r_bar_j=J.q[:3]
    r_bar_i=I.q[:3]
    r_j=np.dot(A_j,r_bar_j)
    r_i=np.dot(A_i,r_bar_i)
    
    r_dot_i=I.r_dot
    r_dot_j=J.r_dot
    
    s_bar_j=J.s_bar
    s_bar_i=I.s_bar    

    B_dot_i=build_B(p_dot_i,a_bar_i)    
    B_dot_j=build_B(p_dot_j,a_bar_j)
    
    B_p_dot_j__s_bar_j=build_B(p_dot_j,s_bar_j)
    B_p_dot_i__s_bar_i=build_B(p_dot_i,s_bar_i)
    
    B_i=build_B(p_i,a_bar_i)
    B_j=build_B(p_j,a_bar_j)
    
    B_p_j__s_bar_j=build_B(p_j,s_bar_j)
    B_p_i__s_bar_i=build_B(p_i,s_bar_i)
    
    a_dot_i=np.dot(B_i,p_dot_i)
    a_dot_j=np.dot(B_j,p_dot_j)
    
    a_i_T=np.transpose(a_i)
    a_j_T=np.transpose(a_j)
    
    a_dot_i_T=np.transpose(a_dot_i)
    
    d_ij_T=np.transpose(r_j+np.dot(A_j,s_bar_j)-r_i-np.dot(A_i,s_bar_i))
    d_ij_dot=r_dot_j+np.dot(B_p_j__s_bar_j,p_dot_j)-r_dot_i-np.dot(B_p_i__s_bar_i,p_dot_i)


    gamma_DP2=np.dot(np.dot(-a_i_T,B_p_dot_j__s_bar_j),p_dot_j)+np.dot(np.dot(a_i_T,B_p_dot_i__s_bar_i),p_dot_i)-np.dot(np.dot(d_ij_T,B_dot_i),p_dot_i)-2*np.dot(a_dot_i_T,d_ij_dot)+ddf

    return float(gamma_DP2)
#%%
def DP2_PHI_partials(X):
      #Set bodys i and j 
    for i in X:
        if i.name == "body i":
            i_index=X.index(i)
        if i.name == "body j":
            j_index=X.index(i) 
            
    I=X[i_index]
    J=X[j_index]  
    #Set a and A for bodies i and j
    a_bar_i=I.a_bar
    a_bar_j=J.a_bar
    A_i=build_A(I.q[3:])
    A_j=build_A(J.q[3:])
    
    #Performe some transposes 
    a_bar_i_T=np.transpose(a_bar_i)
    a_bar_j_T=np.transpose(a_bar_j)
    A_i_T=np.transpose(A_i)
    A_j_T=np.transpose(A_j)

    p_i=I.q[3:]
    p_j=J.q[3:]
    p_dot_i=I.p_dot
    p_dot_j=J.p_dot

    s_bar_j=J.s_bar
    s_bar_i=I.s_bar
    
    r_i=I.q[:3]
    r_j=J.q[:3]
    
    a_j=np.dot(J.A_rotation,a_bar_j)
    a_i=np.dot(I.A_rotation,a_bar_i)
    a_i_T=np.transpose(a_i)


    a_bar_tilde_i=tilde(a_bar_i)

    a_bar_tilde_j=tilde(a_bar_j)
    a_j_T=np.transpose(a_j)
    
    
    B_i=build_B(p_i,a_bar_i)
    B_j=build_B(p_j,s_bar_j)

    
    d_ij=(r_j+np.dot(A_j,s_bar_j)-r_i-np.dot(A_i,s_bar_i))
    d_ij_T=np.transpose(d_ij)    
    
    PHI_DP1_p_i=np.dot(d_ij_T,B_i)
    PHI_DP1_p_j=np.dot(a_i_T,B_j)   
    
    PHI_DP1_r_i=-a_i_T
    PHI_DP1_r_j=a_i_T    
    
    return PHI_DP1_r_i[0], PHI_DP1_r_j[0],PHI_DP1_p_i[0],PHI_DP1_p_j[0]
