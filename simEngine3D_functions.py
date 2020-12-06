import numpy as np
from simEngine3D_dataload import data_file,body
import pandas as pd

from CD import *
from DP1 import *
#%%

# def df(t):
#     return -np.pi/2. * np.sin(2.* t) * np.cos(np.pi/4.* np.cos(2.* t))
#     # return 0

# def ddf(t):
#     return -np.pi/4. * (4.* np.cos(2.* t) * np.cos(np.pi/4.* np.cos(2* t)) + \
#             np.pi * np.sin(2.* t) * np.sin(2.* t) * np.sin(np.pi/4.* np.cos(2.* t)))
#     # return 0
    
#%%
# calculate a_tilde
def tilde(a):
    ax, ay, az = float(a[0]), float(a[1]), float(a[2])

    t = np.array([[0, -az, ay],
                  [az, 0, -ax],
                  [-ay, ax, 0]])
    return t
#%%
def build_M(SYS):
    M=np.zeros((3*SYS.nb,3*SYS.nb))
    for i in range(0,SYS.nb):
        for j in range(0,3):
            M[3*i+j,3*i+j]=SYS.bodies[i].mass
    return M
#%%
def build_JP(SYS):
    J_P=np.zeros((4*SYS.nb,4*SYS.nb))
    for i in range(0,SYS.nb):
        G=build_G(SYS.bodies[i].q[3:])
        J_bar=J_bar_bar(SYS.bodies[i])
        J_P[0+(4*i):4+4*i,0+(4*i):4+4*i]=4*np.dot(np.dot(np.transpose(G),J_bar),G)    

    return J_P

#%%
def build_P(SYS):
    P=np.zeros([SYS.nb,SYS.nb*4])
    for i in range(0,SYS.nb):
        P[i,:]=np.transpose(SYS.bodies[i].q[3:])
        
    return P

#%%
def build_F(SYS):
    F=np.zeros([3*SYS.nb])
    for i in range(0,SYS.nb):
        F[0+3*i:3+3*i]=SYS.bodies[i].gravity*SYS.bodies[i].mass
           
    F.shape=(3,1)
    return F
#%%
def build_tau_hat(SYS):
    tau_hat=np.zeros((4*SYS.nb,1))
    for i in range(0,SYS.nb):
        G_dot=build_G(SYS.bodies[i].p_dot)
        J_bar=J_bar_bar(SYS.bodies[i])        
        tau_hat[0+4*i:4+4*i]=8*np.dot(np.dot(np.dot(G_dot.T,J_bar),G_dot),SYS.bodies[i].q[3:])        
        
    return tau_hat
#%%
def build_gamma_p(SYS):
    gamma_p=np.zeros((SYS.nb,1))
    for i in range(0,SYS.nb):
        gamma_p[i]=-2*np.dot(SYS.bodies[0].p_dot.T,SYS.bodies[0].p_dot) #slide 473
    
    return gamma_p
    
#%%
def J_bar_bar(X):
    m=X.mass
    J_bar=np.zeros((3,3))
    L=X.dimensions[0]*2    
    b=X.dimensions[1]/2
    c=X.dimensions[2]/2
    
    J_bar[0,0]=1/12*m*(b**2+c**2)
    J_bar[1,1]=1/12*m*(L**2+c**2)
    J_bar[2,2]=1/12*m*(L**2+b**2)
    
    return J_bar
    


#%%
def body_count(X):
    count=0
    for i in X:
        if i.ground == 'True':
            count=count
        else:
            count=count+1
    return count
    
#%%
def compute_accel(SYS):
    
    SYS.M=build_M(SYS)
    SYS.J_P=build_JP(SYS)
    SYS.P=build_P(SYS)
    SYS.F=build_F(SYS)
    SYS.tau_hat=build_tau_hat(SYS)
    SYS.gamma_p=build_gamma_p(SYS)
    # SYS.gamma=calc_gamma(SYS.bodies,SYS.constraints,SYS.ddf)
    SYS.gamma=calc_gamma(SYS)    
    
    SYS.gamma_hat=SYS.gamma[0:SYS.nc] #remove euler constraint here
    SYS.gamma_hat.shape=(SYS.nc,1)
    
    
    SYS.phi_q=calc_partials(SYS.bodies,SYS.constraints)    
    phi_r=SYS.phi_q[:,0:3*SYS.nb]
    phi_p=SYS.phi_q[:,SYS.nb*3:]
    
    phi_r=phi_r[0:SYS.nc] #don't consider euler constraint
    phi_p=phi_p[0:SYS.nc]
    
    
    RHS=np.zeros([3*SYS.nb+4*SYS.nb+SYS.nb+SYS.nc,1])
    RHS[0:3*SYS.nb]=SYS.F
    RHS[3*SYS.nb:3*SYS.nb+4*SYS.nb]=SYS.tau_hat
    RHS[3*SYS.nb+4*SYS.nb:3*SYS.nb+4*SYS.nb+SYS.nb]=SYS.gamma_p
    RHS[3*SYS.nb+4*SYS.nb+SYS.nb:3*SYS.nb+4*SYS.nb+SYS.nb+SYS.nc]=SYS.gamma_hat

    #sizes taken from slide 667 of all fall 2020 lectures       
    LHS=np.zeros([3*SYS.nb+4*SYS.nb+SYS.nb+SYS.nc,3*SYS.nb+4*SYS.nb+SYS.nb+SYS.nc])
    LHS[0:3*SYS.nb,0:3*SYS.nb]=SYS.M
    LHS[0:3*SYS.nb,3*SYS.nb+4*SYS.nb+SYS.nb:3*SYS.nb+4*SYS.nb+SYS.nb+SYS.nc]=phi_r.T       
    LHS[3*SYS.nb:3*SYS.nb+4*SYS.nb,3*SYS.nb:3*SYS.nb+4*SYS.nb]=SYS.J_P    
    LHS[3*SYS.nb:3*SYS.nb+4*SYS.nb,3*SYS.nb+4*SYS.nb:3*SYS.nb+4*SYS.nb+SYS.nb]=SYS.P.T
    LHS[3*SYS.nb:3*SYS.nb+4*SYS.nb,3*SYS.nb+4*SYS.nb+SYS.nb:3*SYS.nb+4*SYS.nb+SYS.nb+SYS.nc]=phi_p.T
    LHS[3*SYS.nb+4*SYS.nb:3*SYS.nb+4*SYS.nb+SYS.nb,3*SYS.nb:3*SYS.nb+4*SYS.nb]=SYS.P
    LHS[3*SYS.nb+4*SYS.nb+SYS.nb:3*SYS.nb+4*SYS.nb+SYS.nb+SYS.nc,0:3*SYS.nb]=phi_r
    LHS[3*SYS.nb+4*SYS.nb+SYS.nb:3*SYS.nb+4*SYS.nb+SYS.nb+SYS.nc,3*SYS.nb:3*SYS.nb+4*SYS.nb]=phi_p       
    
    unknowns=np.dot(np.linalg.inv(LHS),RHS)    
    
    ddr=unknowns[0:3*SYS.nb]
    ddp=unknowns[3*SYS.nb:3*SYS.nb+4*SYS.nb]    
    lambda_p=unknowns[3*SYS.nb+4*SYS.nb:3*SYS.nb+4*SYS.nb+SYS.nb]
    lagrange=unknowns[3*SYS.nb+4*SYS.nb+SYS.nb:3*SYS.nb+4*SYS.nb+SYS.nb+SYS.nc]
    
    columns_to_update=["ddr0","ddp0","lambda_p0","lagrange0"]    
    values=[ddr,ddp,lambda_p,lagrange]
    
    for i in range(0,len(columns_to_update)):
        col=columns_to_update[i]
        val=values[i]
        s_n=SYS.outputs[col].copy()
        s_n[SYS.n]=val
        SYS.outputs[col]=s_n
    
    # for i,j in update_dict.items():
    #     col=i    
    #     values_slice=j
    #     s_n=SYS.outputs[col].copy()
    #     s_n[SYS.n]=unknowns[values_slice]       
        
        # SYS.outputs[col]=s_n          
    
    # columns_to_update=["ddr0","ddp0","lambda_p0","lagrange0"]
    # update_values=[slice(0,3*SYS.nb,1),slice(3*SYS.nb,3*SYS.nb+4*SYS.nb,1),slice(3*SYS.nb+4*SYS.nb,3*SYS.nb+4*SYS.nb+SYS.nb,1),slice(3*SYS.nb+4*SYS.nb+SYS.nb,3*SYS.nb+4*SYS.nb+SYS.nb+SYS.nc,1)]
    # for i in range(0,len(columns_to_update)):
    #     col=columns_to_update[i]
    #     values_slice=update_values[i]
    #     s_n=SYS.outputs[col].copy()
    #     s_n[SYS.n]=unknowns[values_slice]       
        
    #     SYS.outputs[col]=s_n        
        

    # SYS.outputs.loc[SYS.n,"ddr0"]=unknowns[0:3*SYS.nb]
    # SYS.outputs.loc[SYS.n,"ddp0"]=unknowns[3*SYS.nb:3*SYS.nb+4*SYS.nb]  
    # SYS.outputs.loc[SYS.n,"lambda_p0"]= unknowns[3*SYS.nb+4*SYS.nb:3*SYS.nb+4*SYS.nb+SYS.nb]
    # SYS.outputs.loc[SYS.n,"lagrange0"]=unknowns[3*SYS.nb+4*SYS.nb+SYS.nb:3*SYS.nb+4*SYS.nb+SYS.nb+SYS.nc]
    
#%%

def build_DF(SYS,time_list):
    outputs=pd.DataFrame(time_list,columns=["time"])
    SYS.cols=["r","dr","ddr","p","dp","ddp","lambda_p","lagrange"]
    for j in range(0,SYS.nb+1):
        for i in SYS.cols:
            if "r" in i:
                outputs[i+str(j)]=[np.zeros([3,1])]*len(outputs)
            if "p" in i:
                outputs[i+str(j)]=[np.zeros([4,1])]*len(outputs)
            if "lambda_p" in i:
                outputs[i+str(j)]=[np.zeros([SYS.nb,1])]*len(outputs)
            if "lagrange" in i:
                outputs[i+str(j)]=[np.zeros([SYS.nc,1])]*len(outputs)
            
        
    return outputs

    
#%%
def update_level0(SYS):
    update_cols=[]
    for i in range(0,SYS.nb+1):
        update_i=[]
        for j in SYS.cols:
            if len(j) ==1:
                update_i.append(j+str(i))
        update_cols.append(update_i)

    for i in update_cols:
        body=SYS.bodies[int(i[0][-1])]
        for j in i:
            if "r" in j:
                r0_n=SYS.outputs[j].copy()
                r0_n[SYS.n]=body.q[:3].copy()
                SYS.outputs[j]=r0_n
                # SYS.outputs.loc[SYS.n,j]=body.q[:3]
            else:
                p0_n=SYS.outputs[j].copy()
                p0_n[SYS.n]=body.q[3:].copy()
                SYS.outputs[j]=p0_n                
                # SYS.outputs.loc[SYS.n,j]=body.q[3:]

#%%
def update_level1(SYS):
    update_cols=[]
    for i in range(0,SYS.nb+1):
        update_i=[]
        for j in SYS.cols:
            if len(j) ==2:
                update_i.append(j+str(i))
        update_cols.append(update_i)

    for i in update_cols:
        body=SYS.bodies[int(i[0][-1])]
        for j in i:
            if "r" in j:
                r_dot=body.r_dot
                r_dot.shape=(3,1)
                
                s_n=SYS.outputs[j].copy() #foolishly updated pandas so I had to change indexing method
                s_n[SYS.n]=r_dot      
                SYS.outputs[j]=s_n
                
                # SYS.outputs.loc[SYS.n,j]=r_dot
            else:
                p_dot=body.p_dot
                p_dot.shape=(4,1)      
                
                s_n=SYS.outputs[j].copy()
                s_n[SYS.n]=p_dot      
                SYS.outputs[j]=s_n
                # SYS.outputs.loc[SYS.n,j]=body.p_dot   
                
#%%
def update_level2(SYS):
    update_cols=[]
    for i in range(0,SYS.nb+1):
        update_i=[]
        for j in SYS.cols:
            if len(j) ==3:
                update_i.append(j+str(i))
        update_cols.append(update_i)

    for i in update_cols:
        body=SYS.bodies[int(i[0][-1])]
        for j in i:
            if "r" in j:
                r_ddot=body.r_ddot.copy()
                r_ddot.shape=(3,1)
                
                s_n=SYS.outputs[j].copy() #foolishly updated pandas so I had to change indexing method
                s_n[SYS.n]=r_ddot      
                SYS.outputs[j]=s_n

            else:
                p_ddot=body.p_ddot.copy()
                p_ddot.shape=(4,1)          
                s_n=SYS.outputs[j].copy() #foolishly updated pandas so I had to change indexing method
                s_n[SYS.n]=p_ddot    
                SYS.outputs[j]=s_n
                
                             
#%%
def BDF(SYS,n):
    
    tol=1e-3
    max_iter=50
    h=SYS.h
    time=SYS.outputs.time[n]

    if n < 2: #use BDF order 1
    
        beta_0 = 1
        alpha1 = -1
    
        sol_1=SYS.outputs.iloc[n-1]        

        C_r_dot = -alpha1*np.array(sol_1.dr0)
        C_r_dot.shape=(3,1)
        C_r = -alpha1*np.array(sol_1.r0) + beta_0*h*C_r_dot
        C_p_dot = -alpha1*np.array(sol_1.dp0)
        C_p_dot.shape=(4,1)
        C_p = -alpha1*np.array(sol_1.p0) + beta_0*h*C_p_dot
    else: #use BDF order 2
        beta_0 = 2/3
        alpha1 = -4/3
        alpha2 = 1/3
    
        sol_1=SYS.outputs.iloc[n-1]               
        sol_2=SYS.outputs.iloc[n-2]     
        
        sol_2.dp0.shape=(4,1)
                  
        C_r_dot = -alpha1*np.array(sol_1.dr0) - alpha2*np.array(sol_2.dr0)
        C_r=-alpha1*np.array(sol_1.r0) - alpha2*np.array(sol_2.r0) + beta_0*h*C_r_dot
        C_p_dot=-alpha1*np.array(sol_1.dp0) - alpha2*np.array(sol_2.dp0)
        C_p_dot.shape=(4,1)        
        C_p = -alpha1*np.array(sol_1.p0) - alpha2*np.array(sol_2.p0) + beta_0*h*C_p_dot 
    
    

    ddr=np.array(sol_1.ddr0)
    ddp=np.array(sol_1.ddp0)
    lambda_p=np.array(sol_1.lambda_p0)
    lagrange=np.array(sol_1.lagrange0)
    
    ddr.shape=(3,1)
    ddp.shape=(4,1)
    
    count=0
    error=1
    while error > tol and count < max_iter:
        r=C_r + beta_0**2*h**2*ddr
        dr=C_r_dot + beta_0*h*ddr
        p=C_p + beta_0**2*h**2*ddp
        dp=C_p_dot + beta_0*h*ddp
        
        update_bodies(SYS,r,p,dr,dp,ddr,ddp)

        M=build_M(SYS)
        J_p = build_JP(SYS)
        P = build_P(SYS)
        F = build_F(SYS)

        tau_hat=build_tau_hat(SYS)

        
                
        SYS.phi_q=calc_partials(SYS.bodies,SYS.constraints)    
        phi_r=SYS.phi_q[:,0:3*SYS.nb]
        phi_p=SYS.phi_q[:,SYS.nb*3:]
        
        phi_r=phi_r[0:SYS.nc] #don't consider euler constraint
        phi_p=phi_p[0:SYS.nc]    
        
        PHI=calc_phi(SYS.bodies,SYS.constraints,time)
        PHI=np.array(PHI[0:SYS.nc]) #remove euler parameterization constraint        
        PHI.shape=(SYS.nc,1)
        
        PHI_P=1/2*(p.T @ p) - 1/2           
                  
        ddr.shape=(3*SYS.nb,1)
        ddp.shape=(4*SYS.nb,1)
        lambda_p.shape=(SYS.nb,1)
        lagrange.shape=(SYS.nc,1)
        g=np.zeros([8*SYS.nb+SYS.nc,1])
        g[0:3*SYS.nb]=M @ ddr + phi_r.T @ lagrange - F
        g[3*SYS.nb:3*SYS.nb+4*SYS.nb] =J_p @ ddp + phi_p.T @ lagrange + P.T @ lambda_p - tau_hat
        g[3*SYS.nb+4*SYS.nb:3*SYS.nb+4*SYS.nb+SYS.nb]=PHI_P/(beta_0**2*h**2)
        g[3*SYS.nb+4*SYS.nb+SYS.nb:3*SYS.nb+4*SYS.nb+SYS.nb+SYS.nc]=PHI/(beta_0**2*h**2)
        
        
        if count == 0:
            PSI_1=np.concatenate((M,np.zeros((3*SYS.nb,4*SYS.nb)),np.zeros((3*SYS.nb,SYS.nb)),phi_r.T),axis=1)
            PSI_2=np.concatenate((np.zeros((4*SYS.nb,3*SYS.nb)),J_p,P.T,phi_p.T),axis=1)
            PSI_3=np.concatenate((np.zeros((SYS.nb,3*SYS.nb)),P,np.zeros((SYS.nb,SYS.nb)),np.zeros((SYS.nb,SYS.nc))),axis=1)
            PSI_4=np.concatenate((phi_r,phi_p,np.zeros((SYS.nc,SYS.nb)),np.zeros((SYS.nc,SYS.nc))),axis=1)            
                
            PSI=np.concatenate((PSI_1,PSI_2,PSI_3,PSI_4),axis=0)
            
            

        delta_z=-np.matmul(np.linalg.inv(PSI),g) 
        
        ddr=ddr+delta_z[0:3*SYS.nb]
        ddp=ddp+delta_z[3*SYS.nb:3*SYS.nb+4*SYS.nb]
        lambda_p=lambda_p+delta_z[3*SYS.nb+4*SYS.nb:3*SYS.nb+4*SYS.nb+SYS.nb]
        lagrange=lagrange+delta_z[3*SYS.nb+4*SYS.nb+SYS.nb:3*SYS.nb+4*SYS.nb+SYS.nb+SYS.nc]
        

        
                                  
        error=np.linalg.norm(delta_z)                      
        count=count+1   
    
    update_bodies(SYS,r,p,dr,dp,ddr,ddp) #update object with final values  
    update_level0(SYS) #update position in output object
    update_level1(SYS) #update velocity in output object
    update_level2(SYS) #update velocity in output object    
    
            
        
        #%%
def update_bodies(SYS,r,p,dr,dp,ddr,ddp):
    r.shape=(3*SYS.nb,1)
    dr.shape=(3*SYS.nb,1)
    p.shape=(4*SYS.nb,1)
    dp.shape=(4*SYS.nb,1)        
    ddr.shape=(3*SYS.nb,1)
    ddp.shape=(4*SYS.nb,1)
    
    for i in range(0,SYS.nb):
        SYS.bodies[i].q[:3]=r[3*i:3*i+3]
        SYS.bodies[i].q[3:]=p[4*i:4*i+4]
        SYS.bodies[i].r_dot=dr[3*i:3*i+3]
        SYS.bodies[i].p_dot=dp[4*i:4*i+4]
        SYS.bodies[i].r_ddot[:3]=ddr[3*i:3*i+3]        
        SYS.bodies[i].p_ddot=ddp[4*i:4*i+4]              

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
#%%
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
    

#%%
def check_euler(X,cl,time,tol):
    PHI_P=1/2*(X[1].q[3:].T @ X[1].q[3:]) - 1/2
    if PHI_P > tol:
        print("initial conditions do not satisfy PHI_P=0")
    else:
        print("initial conditions satisfy PHI_P=0")   
    return

#%%

def calc_phi(X,cl,t):
    phi_values=[]
    for i in cl:
        if i.type.strip(" ' ") == 'CD':
            phi_values.append(CD_phi(X,i,eval(i.f)))
        elif i.type.strip(" ' ") == 'DP1':
            phi_values.append(DP1_phi(X,i,eval(i.f)))
            
    body_count=0
    phi_p=[]
    for k in X:
        if k.ground=='False':
            phi_p.append(0.5*(np.dot(k.q[3:].T[0],k.q[3:].T[0])-1))
            body_count=body_count+1
    
    for i in range(0,len(phi_p)):
        phi_values.append(phi_p[i])
        
        
    phi=np.zeros(len(phi_values))
    for i in range(0,len(phi_values)):
        phi[i]=phi_values[i]

    return phi



#%%
def calc_partials(X,cl):
    phi_partials_values=[]
    for i in cl:
            if i.type.strip(" ' ") == 'CD':
                CD_i=CD_PHI_partials(X,i)
                phi_partials_values.append(CD_i)
            if i.type.strip(" ' ") == 'DP1':
                DP_1_i=DP1_PHI_partials(X,i)
                phi_partials_values.append(DP_1_i)    
        
    body_count=0
    for k in X:
        if k.name=='body i':
            p_i=k.q[3:]
        if k.ground=='False':
            body_count=body_count+1
        
    phi_q=np.zeros([len(cl)+body_count,7*body_count])        
    for i in range(0,len(phi_partials_values)):
        phi_q[i,:]=phi_partials_values[i]
        
    p_i.shape=(4,1)
    phi_p=np.concatenate((np.zeros([3,1]),p_i))
    phi_q[-1,:]=phi_p.T
        
    return phi_q

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

def calc_nue(SYS):
    X=SYS.bodies
    cl=SYS.constraints
    t=SYS.time
    
    nue_values=[]
    for i in cl:
            if i.type.strip(" ' ") == 'CD':
                nue_values.append(CD_nue(X,i,SYS.df(t)))
            if i.type.strip(" ' ") == 'DP1':
                nue_values.append(DP1_nue(X,i,SYS.df(t)))            

    nu_p=[]            
    body_count=0
    for k in X:
        if k.ground=='False':
            nu_p.append(0)
            body_count=body_count+1
            
    for i in range(0,len(nu_p)):
        nue_values.append(nu_p[i])
        
    nue=np.zeros(len(nue_values))
    for i in range(0,len(nue_values)):
        nue[i]=nue_values[i]
    
    return nue

#%%

def calc_gamma(SYS):
    X=SYS.bodies
    cl=SYS.constraints
    t=SYS.time
    
    gamma_values=[]
    for i in cl:
            if i.type.strip(" ' ") == 'CD':
                gamma_values.append(gamma_CD(X,i,SYS.ddf(t)))
            if i.type.strip(" ' ") == 'DP1':
                gamma_values.append(gamma_DP1(X,i,SYS.ddf(t)))
                
    gamma_p=[]            
    body_count=0
    for k in X:
        if k.ground=='False':
            p_dot=k.p_dot
            p_dot.shape=(4)
            gamma_p.append(-np.dot(p_dot,p_dot))
            body_count=body_count+1
            
    for i in range(0,len(gamma_p)):
        gamma_values.append(gamma_p[i])
        
    gamma=np.zeros(len(gamma_values))
    for i in range(0,len(gamma_values)):
        gamma[i]=gamma_values[i]
                    
            
    return gamma

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
#%%
def getJp(p, m, w, L):
    # here L is the length of pendulum
    # in L-RF
    J_bar = 1./12. * m * np.array([2.* w**2, L**2 + w**2, L**2 + w**2]) * np.eye(3)
    G = getG(p)

    return 4. * np.matmul(np.transpose(G), np.matmul(J_bar, G))

#%%
def getG(p):
    e0, e1, e2, e3 = p[0], p[1], p[2], p[3]
    return np.array([[-e1, e0, e3, -e2],
                     [-e2, -e3, e0, e1],
                     [-e3, e2, -e1, e0]])
    
#%%
def getE(p):
    e0, e1, e2, e3 = p[0], p[1], p[2], p[3]
    return np.array([[-e1, e0, -e3, e2],
                     [-e2, e3, e0, -e1],
                     [-e3, -e2, e1, e0]])