import numpy as np
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
def phi_d(X,C,ft):
    for i in X:
        if i.name == "body i":
            i_index=X.index(i)
        if i.name == "body j":
            j_index=X.index(i) 
               
    I=X[i_index]
    J=X[j_index]    
    
    a_bar_i=I.a_bar

    A_i=I.A_rotation
    A_j=J.A_rotation

    r_j=J.q[0:3]
    s_bar_j=J.s_bar
    
    r_i=I.q[0:3]
    s_bar_i=I.s_bar
    

    a_bar_i_T=np.transpose(a_bar_i)
    d_ij=(r_j+np.dot(A_j,s_bar_j)-r_i-np.dot(A_i,s_bar_i))    
    
    phi_d=np.dot(np.transpose(d_ij),d_ij)-ft
             

    return float(phi_d)
