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
def gamma_DP1(X,C,ddf):
        #Set bodys i and j 
    for i in X:
        if i.name == "body i":
            i_index=X.index(i)
        if i.name == "body j":
            j_index=X.index(i) 
            
    I=X[i_index]
    J=X[j_index]  
    
    #Set a and A for bodies i and j
    a_bar_i=C.a_bar_i
    a_bar_j=C.a_bar_j
    A_i=I.A_rotation
    A_j=J.A_rotation
    
    p_i=I.q[3:]
    p_j=J.q[3:]
    p_dot_i=I.p_dot
    p_dot_j=J.p_dot

    a_i=np.dot(A_i,a_bar_i)
    a_j=np.dot(A_j,a_bar_j)

    B_dot_i=build_B(p_dot_i,a_bar_i)    
    B_dot_j=build_B(p_dot_j,a_bar_j)
    
    B_i=build_B(p_i,a_bar_i)
    B_j=build_B(p_j,a_bar_j)
    
    a_dot_i=np.dot(B_i,p_dot_i)
    a_dot_j=np.dot(B_j,p_dot_j)
    
    a_i_T=np.transpose(a_i)
    a_j_T=np.transpose(a_j)
    
    a_dot_i_T=np.transpose(a_dot_i)
    
    if C.name == "driving":
        ddf=ddf
    else:
        ddf=0
    
    gamma_DP1=np.dot(np.dot(-a_i_T,B_dot_j),p_dot_j)-np.dot(np.dot(a_j_T,B_dot_i),p_dot_i)-2*np.dot(a_dot_i_T,a_dot_j)+ddf
    return float(gamma_DP1)

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
    A_i=I.A_rotation
    A_j=J.A_rotation
    
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
def gamma_CD(X,C,ddf):
        #Set bodys i and j 
    for i in X:
        if i.name == "body i":
            i_index=X.index(i)
        if i.name == "body j":
            j_index=X.index(i) 
            
    I=X[i_index]
    J=X[j_index]  
    
    p_dot_i=I.p_dot
    p_dot_j=J.p_dot
    
    s_bar_i=C.s_bar_i
    s_bar_j=C.s_bar_j
    
    c=C.c
    
    c_T=np.transpose(c)

    if C.name == "driving":
        ddf=ddf
    else:
        ddf=0
        

    B_dot_i=build_B(p_dot_i,s_bar_i)
    
    if s_bar_j.all() == 0:
        B_dot_j=np.zeros([3,4])
    else:
        B_dot_j=build_B(p_dot_j,s_bar_j)
    
    gamma_CD=np.dot(np.dot(c_T,B_dot_i),p_dot_i)-np.dot(np.dot(c_T,B_dot_j),p_dot_j)+ddf
    
    return float(gamma_CD)
#%%
def gamma_D(X,ddft):
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
    
    #Performe some transposes 
    a_bar_i_T=np.transpose(a_bar_i)
    a_bar_j_T=np.transpose(a_bar_j)
    A_i_T=np.transpose(A_i)
    A_j_T=np.transpose(A_j)
    

    p_dot_i=I.p_dot
    p_dot_j=J.p_dot
    s_bar_i=I.s_bar
    s_bar_j=J.s_bar
    
    r_dot_j=J.r_dot
    r_j=J.q[:3]
    r_i=I.q[:3]
    r_dot_i=I.r_dot
    
    p_i=I.q[3:]
    p_j=J.q[3:]
    

    a_j=np.dot(J.A_rotation,J.a_bar)
    a_i=np.dot(I.A_rotation,I.a_bar)
    a_i_T=np.transpose(a_i)


    a_bar_tilde_i=tilde(a_bar_i)

    a_bar_tilde_j=tilde(a_bar_j)
    a_bar_tilde_j_T=np.transpose(a_bar_tilde_j)
    
    B_dot_i=build_B(p_dot_i,s_bar_i)
    B_dot_j=build_B(p_dot_j,s_bar_j)
    B_p_j__s_bar_j=build_B(p_j,s_bar_j)
    B_p_i__s_bar_i=build_B(p_i,s_bar_i)
    
    d_ij_T=np.transpose(r_j+np.dot(A_j,s_bar_j)-r_i-np.dot(A_i,s_bar_i))
    d_ij_dot=r_dot_j+np.dot(B_p_j__s_bar_j,p_dot_j)-r_dot_i-np.dot(B_p_i__s_bar_i,p_dot_i)
    

    
    gamma_D=-2*np.dot(np.dot(d_ij_T,B_dot_j),p_dot_j)+2*np.dot(np.dot(d_ij_T,B_dot_i),p_dot_i)-2*np.dot(np.transpose(d_ij_dot),d_ij_dot)+ddft
    
    
    return float(gamma_D)