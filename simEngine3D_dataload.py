import numpy as np
import os


#%% Define new class to hold body data
class body():
    """Object that holds information about bodies"""
    def __init__(self):
        self.name = ""
        self.ID = ""    
        self.mass = np.float()
        self.q = np.array([0,0,0,0,0,0,0])
        self.a_bar= np.array([0,0,0],dtype=np.float64)
        self.A_rotation=np.empty([3,3],dtype=np.float64)
        self.r_dot= np.array([0,0,0],dtype=np.float64)
        self.s_bar = np.array([0,0,0],dtype=np.float64)
        self.normalized= bool()
#%%
def array_dot(V):
    for i in range(0,len(V)-1):
        if len(V) >= 2:
            x=np.dot(V[0],V[1])
            V.pop(0)
            if len(V) >1:
                V.pop(0)
                V.insert(0,x)
    
#%%import data from input file
def data_file():

    #File path to input file and read contents of file into memory
    file=os.getcwd()+"\\input.txt"
    
    with open(file,'r') as f:
        contents=f.readlines()
        
    #%% Create empty lists
    
    bodys=[]
    body_list=[]
    body_ob_list=[]
      
        
    #%%
    bodys_index=contents.index('bodies:\n') #Know where to start reading in body information
    constraints_index=contents.index('constraints:\n') #Know where to start reading in constraint information.
    
    #%%Read all lines between bodys and constraint information
    for i in range(bodys_index,constraints_index-1):
        bodys.append(contents[i].strip())
        
    #%%
    start_count_index=0
    end_count_index=0
    count=bodys.count("{") #Count how many left squigly brakets there are, this gives how many bodies there are in the system
    
    #Create list of body objects                 
    for i in range(0,count):
        body_ob_list.append(body())
    
    #loop through all bodies
    for k in range(0,count):
        st=bodys.index("{",start_count_index) #Get first body information
        en=bodys.index("}",end_count_index) #Stop body information once "}" is found
    
        #Create list of each body's information
        x=[]
        for i in range(st+1,en):
            y=bodys[i].partition(":")
            x.append({y[0].strip("\'").strip().strip(" ' "):y[2].strip()})
        
        #assign values to object
        for i in x:
            if  list(i.keys())[0] == "name":
                body_ob_list[k].name= list(i.values())[0].strip(" ' ")
            if  list(i.keys())[0] == "ID":
                body_ob_list[k].ID= int(list(i.values())[0])
            if  list(i.keys())[0] == "mass":
                body_ob_list[k].mass= np.float(list(i.values())[0])
            if  list(i.keys())[0] == "q(x,y,z,e0,e1,e2,e3)":
                a=list(i.values())[0].strip()
                a=a.strip("[")
                a=a.strip("]")
                a=a.split(",")
                q=np.empty([7,1],dtype=float)
                for n in range(0,len(a)):
                    q[n]=float(a[n])
                body_ob_list[k].q= q
            if  list(i.keys())[0] == "a_bar":
                a=list(i.values())[0].strip()
                a=a.strip("[")
                a=a.strip("]")
                a=a.split(",")
                a_=np.empty([3,1],dtype=float)
                for n in range(0,len(a)):
                    a_[n]=float(a[n])
                body_ob_list[k].a_bar= a_
            if  list(i.keys())[0] == 'p_dot':          
                a=list(i.values())[0].strip()
                a=a.strip("[")
                a=a.strip("]")
                a=a.split(",")       
                q=np.empty([4,1],dtype=float)                
                for n in range(0,len(a)):
                    q[n]=float(a[n])
                body_ob_list[k].p_dot = q
            if  list(i.keys())[0] == 's_bar':          
                a=list(i.values())[0].strip()
                a=a.strip("[")
                a=a.strip("]")
                a=a.split(",")       
                q=np.empty([3,1],dtype=float)                
                for n in range(0,len(a)):
                    q[n]=float(a[n])
                body_ob_list[k].s_bar = q            
            if  list(i.keys())[0] == 'r_dot':          
                a=list(i.values())[0].strip()
                a=a.strip("[")
                a=a.strip("]")
                a=a.split(",")       
                q=np.empty([3,1],dtype=float)                
                for n in range(0,len(a)):
                    q[n]=float(a[n])
                body_ob_list[k].r_dot = q        
            if  list(i.keys())[0] == "Need normalized":
                body_ob_list[k].normalized= list(i.values())[0].strip(" ' ")                            
        
        body_list.append(x)
        start_count_index=st+1
        end_count_index=en+1
        
    for i in body_ob_list:
        if i.normalized == 'True':
            i.q[3:]=i.q[3:] / np.linalg.norm(i.q[3:])
            
        
    #Assemble rotation matrix A using the values e0-e3 provided
    A_list=[]
    for i in range(0,len(body_ob_list)):
        A=np.empty([3,3])

        for k in range(0,len(A)):
            A[k,k]=2*((body_ob_list[i].q[3]**2)+(body_ob_list[i].q[4+k]**2)-0.5)
        A[0,1]=2*(body_ob_list[i].q[4]*body_ob_list[i].q[5]-body_ob_list[i].q[3]*body_ob_list[i].q[6])
        A[0,2]=2*(body_ob_list[i].q[4]*body_ob_list[i].q[6]+body_ob_list[i].q[3]*body_ob_list[i].q[5])
        
        A[1,0]=2*((body_ob_list[i].q[4]*body_ob_list[i].q[5]+body_ob_list[i].q[3]*body_ob_list[i].q[6]))
        A[1,2]=2*(body_ob_list[i].q[5]*body_ob_list[i].q[6]-body_ob_list[i].q[3]*body_ob_list[i].q[4])
        
        A[2,0]=2*(body_ob_list[i].q[4]*body_ob_list[i].q[6]-body_ob_list[i].q[3]*body_ob_list[i].q[5])
        A[2,1]=2*((body_ob_list[i].q[5]*body_ob_list[i].q[6]+body_ob_list[i].q[3]*body_ob_list[i].q[4]))
    
        A_list.append(A)
    
    for i in range(0,len(A_list)):
        body_ob_list[i].A_rotation=A_list[i]
        

            
    
    return body_ob_list #returns list of bodies with corisponding information


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
def DP1_PHI_partials(X,C):
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
    
    #Performe some transposes 
    a_bar_i_T=np.transpose(a_bar_i)
    a_bar_j_T=np.transpose(a_bar_j)
    A_i_T=np.transpose(A_i)
    A_j_T=np.transpose(A_j)

    p_i=I.q[3:]
    p_j=J.q[3:]
    p_dot_i=I.p_dot
    p_dot_j=J.p_dot

    a_j=np.dot(J.A_rotation,J.a_bar)
    a_i=np.dot(I.A_rotation,I.a_bar)
    a_i_T=np.transpose(a_i)


    a_bar_tilde_i=tilde(a_bar_i)

    a_bar_tilde_j=tilde(a_bar_j)
    a_j_T=np.transpose(a_j)
    
    
    B_i=build_B(p_i,a_bar_i)
    B_j=build_B(p_j,a_bar_j)
    PHI_DP1_p_i=np.dot(a_j_T,B_i)
    PHI_DP1_p_j=np.dot(a_i_T,B_j)   
    
    PHI_DP1_r_i=0
    PHI_DP1_r_j=0
    
    return PHI_DP1_r_i, PHI_DP1_r_j,PHI_DP1_p_i[0],PHI_DP1_p_j[0]
#%%
def CD_PHI_partials(X,C):
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
    
    #Performe some transposes 
    a_bar_i_T=np.transpose(a_bar_i)
    a_bar_j_T=np.transpose(a_bar_j)
    A_i_T=np.transpose(A_i)
    A_j_T=np.transpose(A_j)
    
    
    #define G and p_dot

    p_dot_i=I.p_dot
    p_dot_j=J.p_dot
    
    p_i=I.q[3:]
    p_j=J.q[3:]
    
    s_bar_i=C.s_bar_i
    s_bar_j=C.s_bar_j
    
    a_j=np.dot(J.A_rotation,J.a_bar)
    a_i=np.dot(I.A_rotation,I.a_bar)
    a_i_T=np.transpose(a_i)


    a_bar_tilde_i=tilde(a_bar_i)

    a_bar_tilde_j=tilde(a_bar_j)
    a_bar_tilde_j_T=np.transpose(a_bar_tilde_j)
    
    B_i=build_B(p_i,s_bar_i)
    if s_bar_j.all()==0:
        B_j=np.zeros([3,4])
    else:
        B_j=build_B(p_j,s_bar_j)
    c=C.c
    c_T=np.transpose(c)
    
    first=-c_T
    second=np.dot(-c_T,B_i)
    third=c_T
    fourth=np.dot(c_T,B_j)
    
    
    return first[0], second[0], third[0], fourth[0]
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
    A_i=I.A_rotation
    A_j=J.A_rotation
    
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
    
    a_j=np.dot(J.A_rotation,J.a_bar)
    a_i=np.dot(I.A_rotation,I.a_bar)
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

#%%
def D_PHI_partials(X):
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

    p_i=I.q[3:]
    p_j=J.q[3:]
    p_dot_i=I.p_dot
    p_dot_j=J.p_dot

    s_bar_j=J.s_bar
    s_bar_i=I.s_bar
    
    r_i=I.q[:3]
    r_j=J.q[:3]
    
    a_j=np.dot(J.A_rotation,J.a_bar)
    a_i=np.dot(I.A_rotation,I.a_bar)
    a_i_T=np.transpose(a_i)


    a_bar_tilde_i=tilde(a_bar_i)
    a_bar_tilde_j=tilde(a_bar_j)
    a_j_T=np.transpose(a_j)
    

    B_p_i__s_bar_i=build_B(p_i,s_bar_i)
    B_p_j__s_bar_j=build_B(p_j,s_bar_j)
    
    d_ij=(r_j+np.dot(A_j,s_bar_j)-r_i-np.dot(A_i,s_bar_i))
    d_ij_T=np.transpose(d_ij)    
    
    PHI_DP1_p_i=-2*np.dot(d_ij_T,B_p_i__s_bar_i)
    PHI_DP1_p_j=2*np.dot(d_ij_T,B_p_j__s_bar_j)
    PHI_DP1_r_i=-2*d_ij_T
    PHI_DP1_r_j=2*d_ij_T    
    
    return PHI_DP1_r_i[0], PHI_DP1_r_j[0],PHI_DP1_p_i[0],PHI_DP1_p_j[0]