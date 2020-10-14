import os
import numpy as np

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


file="C:\\Users\\Logan\\Desktop\\simEngine3D\\input.txt"

def constraints_in():
    with open(file,'r') as f:
        contents=f.readlines()
        
        
    constraints_index=contents.index('constraints:\n') #Know where to start reading in body information
    end_index=contents.index(']\n',constraints_index) #Know where to start reading in body information
    
    constraints=[]    
    start_count_index=0
    end_count_index=0                      
                              
    
    for i in range(constraints_index,end_index):
        constraints.append(contents[i].strip())    
        
        
    count=constraints.count("{") #Count how many left squigly brakets there are, this gives how many bodies there are in the system    

    constraint_list=[]
    for i in range(0,count):
        constraint_list.append(constraint())
        
    C=[]
    for k in range(0,count):
        st=constraints.index("{",start_count_index) #Get first body information
        en=constraints.index("}",end_count_index) #Stop body information once "}" is found

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
                                                        
        start_count_index=st+1
        end_count_index=en+1                                       

          
          
    return constraint_list