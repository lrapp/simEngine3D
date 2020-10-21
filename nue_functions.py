#nue functions

import numpy as np

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
def DP1_nue(X,C,dft):

    if C.name == "driving":
        dft=dft
    else:
        dft=0

    return dft
#%%
def CD_nue(X,C,dft):

    if C.name == "driving":
        dft=dft
    else:
        dft=0
        
    return dft

#%%
def DP2_nue(X,dft):
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
    A_i=I.A_rotation
    A_j=J.A_rotation


    return dft
#%%
def D_nue(X,dft):
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
    A_i=I.A_rotation
    A_j=J.A_rotation
    


    return dft