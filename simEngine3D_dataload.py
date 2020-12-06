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
        self.dimensions= np.array([0,0,0],dtype=np.float64)        
        self.a_bar= np.array([0,0,0],dtype=np.float64)
        self.r_dot= np.array([0,0,0],dtype=np.float64)
        self.r_ddot= np.zeros([3,1])
        self.s_bar = np.array([0,0,0],dtype=np.float64)
        self.normalized= bool()
        self.p_dot= np.zeros([4,1])
        self.p_ddot= np.zeros([4,1])    
        self.ground=""
        self.type=""

#%%import data from input file
def data_file(file):

    #File path to input file and read contents of file into memory
#    file=os.getcwd()+"\\input.txt"
#    file=os.getcwd()+"\\revJoint.mdl"    
    
    with open(file,'r') as f:
        contents=f.readlines()
        
    #%Create empty lists
    bodys=[]
    body_list=[]
    body_ob_list=[]
      
        
    bodys_index=contents.index('bodies:\n') #Know where to start reading in body information
    constraints_index=contents.index('constraints:\n') #Know where to start reading in constraint information.
    
    #Read all lines between bodys and constraint information
    for i in range(bodys_index,constraints_index-1):
        bodys.append(contents[i].strip())
        
    #
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
                
            if  list(i.keys())[0] == "ground":
                body_ob_list[k].ground= list(i.values())[0].strip(" ' ")

            if  list(i.keys())[0] == "type":
                body_ob_list[k].type= list(i.values())[0].strip(" ' ")                
                
            if  list(i.keys())[0] == "dimensions":
                a=list(i.values())[0].strip()
                a=a.strip("[")
                a=a.strip("]")
                a=a.split(",")
                q=np.zeros((3),dtype=float)
                for n in range(0,len(a)):
                    q[n]=float(a[n])
                body_ob_list[k].dimensions= q  
            if  list(i.keys())[0] == "gravity":
                a=list(i.values())[0].strip()
                a=a.strip("[")
                a=a.strip("]")
                a=a.split(",")
                q=np.zeros((3),dtype=float)
                for n in range(0,len(a)):
                    q[n]=float(a[n])
                body_ob_list[k].gravity= q                                
        
        body_list.append(x)
        start_count_index=st+1
        end_count_index=en+1
        
    for i in body_ob_list:
        if i.normalized == 'True':
            i.q[3:]=i.q[3:] / np.linalg.norm(i.q[3:])
                  
    
    return body_ob_list #returns list of bodies
#%%
class constraint():
    """Object that holds information about bodies"""
    def __init__(self):
        self.name = ""
        self.ID = ""    
        self.type = ""
        self.a_bar_i= np.array([0,0,0],dtype=np.float64)
        self.a_bar_j= np.array([0,0,0],dtype=np.float64)
        self.s_bar_i = np.array([0,0,0],dtype=np.float64)
        self.s_bar_j = np.array([0,0,0],dtype=np.float64)
        self.f=float()
        self.body_i_name=""
        self.body_j_name=""



def constraints_in(file):
    with open(file,'r') as f:
        contents=f.readlines()
        
    constraints_index=contents.index('constraints:\n') #Know where to start reading in constraint information
    end_index=contents.index(']\n',constraints_index) #Know where to stop reading in constraint information
    
    constraints=[]    
    start_count_index=0
    end_count_index=0                      
                              
    
    for i in range(constraints_index,end_index):
        constraints.append(contents[i].strip())    
        
        
    count=constraints.count("{") #Count how many left squigly brakets there are, this gives how many constraints there are in the system    

    constraint_list=[]
    for i in range(0,count):
        constraint_list.append(constraint())
        
    C=[]
    for k in range(0,count):
        st=constraints.index("{",start_count_index) #Get first constraint information
        en=constraints.index("}",end_count_index) #Stop constraint information once "}" is found

        y=[]
        x=[]
        for i in range(st+1,en):
        #Create list of each body's information
            x=[]
            for i in range(st+1,en):
                y=constraints[i].partition(":")
                x.append({y[0].strip("\'").strip().strip(" ' "):y[2].strip()})            
            
                
        for i in x:           
                if  list(i.keys())[0] == "name":
                    constraint_list[k].name= list(i.values())[0].strip(" ' ")            
                if  list(i.keys())[0] == "ID":                    
                    constraint_list[k].ID= int(list(i.values())[0])                    
                if  list(i.keys())[0] == "type":
                    constraint_list[k].type= list(i.values())[0]                        
                if  list(i.keys())[0] == 's_bar_i':          
                    a=list(i.values())[0].strip()
                    a=a.strip("[")
                    a=a.strip("]")
                    a=a.split(",")       
                    q=np.empty([3,1],dtype=float)                
                    for n in range(0,len(a)):
                        q[n]=float(a[n])
                        constraint_list[k].s_bar_i=q
                if  list(i.keys())[0] == 's_bar_j':          
                    a=list(i.values())[0].strip()
                    a=a.strip("[")
                    a=a.strip("]")
                    a=a.split(",")       
                    q=np.empty([3,1],dtype=float)                
                    for n in range(0,len(a)):
                        q[n]=float(a[n])
                        constraint_list[k].s_bar_j=q            
                if  list(i.keys())[0] == 'a_bar_i':          
                    a=list(i.values())[0].strip()
                    a=a.strip("[")
                    a=a.strip("]")
                    a=a.split(",")       
                    q=np.empty([3,1],dtype=float)                
                    for n in range(0,len(a)):
                        q[n]=float(a[n])
                        constraint_list[k].a_bar_i=q                
                if  list(i.keys())[0] == 'a_bar_j':          
                    a=list(i.values())[0].strip()
                    a=a.strip("[")
                    a=a.strip("]")
                    a=a.split(",")       
                    q=np.empty([3,1],dtype=float)                
                    for n in range(0,len(a)):
                        q[n]=float(a[n])
                        constraint_list[k].a_bar_j=q                                        
                if  list(i.keys())[0] == 'c':          
                    a=list(i.values())[0].strip()
                    a=a.strip("[")
                    a=a.strip("]")
                    a=a.split(",")       
                    q=np.empty([3,1],dtype=float)                
                    for n in range(0,len(a)):
                        q[n]=float(a[n])
                        constraint_list[k].c=q           
                if  list(i.keys())[0] == 'f':
                    if list(i.values())[0].strip()==0:
                        constraint_list[k].f=0
                    else:
                        constraint_list[k].f=list(i.values())[0].strip()
                if  list(i.keys())[0] == 'between':
                    bodies=list(i.values())[0].strip().split('and')
                    for ii in bodies:
                        bodies[bodies.index(ii)]=ii.strip(" ' ")
                    constraint_list[k].body_i_name=bodies[0]
                    constraint_list[k].body_j_name=bodies[1]
                        
                                                        
        start_count_index=st+1
        end_count_index=en+1                                       

          
          
    return constraint_list