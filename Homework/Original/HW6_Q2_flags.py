import os
import numpy as np

class constraint():
    """Object that holds information about bodies"""
    def __init__(self):
        self.name = ""
        self.ID = ""    
        self.type = ""
        self.a_bar= np.array([0,0,0],dtype=np.float64)
        self.s_bar = np.array([0,0,0],dtype=np.float64)

#File path to input file and read contents of file into memory
#file=os.getcwd()+"\\input.txt"
file="C:\\Users\\Logan\\Desktop\\ME-751\\HW6\\Q2\\input.txt"
def flags_in():
    with open(file,'r') as f:
        contents=f.readlines()
        
        
    flags_index=contents.index('flags:\n') #Know where to start reading in body information
    bodys_index=contents.index('bodies:\n') #Know where to start reading in body information
    
    flags=[]                          
                              
    for i in range(flags_index,bodys_index-1):
        flags.append(contents[i].strip())    
    
    st=flags.index("{") #Get first body information
    en=flags.index("}") #Stop body information once "}" is found                           
        
    y=[]
    for i in range(st+1,en):
        y=flags[i].split(",")
    return y

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
            
            
            if  list(i.keys())[0] == 's_bar':          
                a=list(i.values())[0].strip()
                a=a.strip("[")
                a=a.strip("]")
                a=a.split(",")       
                q=np.empty([3,1],dtype=float)                
                for n in range(0,len(a)):
                    q[n]=float(a[n])

            for k in y:
              x.append(k.strip(" "))  
        C.append(x)
          
          
    return x

def points_in():
    with open(file,'r') as f:
        contents=f.readlines()
        
        
    points_index=contents.index('points:\n') #Know where to start reading in body information

    
    points=[]   
                              
    
    for i in range(points_index,len(contents)-1):
        points.append(contents[i].strip())    
        
    start_count_index=0
    end_count_index=0
    count=points.count("{") #Count how many left squigly brakets there are, this gives how many bodies there are in the system
    
        
    #loop through all points
    for k in range(0,count):
        st=points.index("{",start_count_index) #Get first body information
        en=points.index("}",end_count_index) #Stop body information once "}" is found
    
        #Create list of each body's information
        x=[]
        for i in range(st+1,en):
            y=points[i].partition(":")
            x.append({y[0].strip("\'").strip().strip(" ' "):y[2].strip()})        
    x=[]
    for i in range(st+1,en):
        y=points[i].partition(":")
        x.append({y[0].strip("\'").strip().strip(" ' "):y[2].strip()})        

    for i in x:
        if  list(i.keys())[0] == 'c':       
            a=list(i.values())[0].strip()
            a=a.strip("[")
            a=a.strip("]")
            a=a.split(",")       
            q=np.empty([3,1],dtype=float)                
            for n in range(0,len(a)):
                q[n]=float(a[n])
                
    return q

def f_in():
    with open(file,'r') as f:
        contents=f.readlines()
        
    f_index=contents.index('ft:\n') #Know where to start reading in body information

    
    f=[]   
                              
    
    for i in range(f_index,contents.index("\t}\n",f_index)+1):
        f.append(contents[i].strip())    
        
    start_count_index=0
    end_count_index=0
    count=f.count("{") #Count how many left squigly brakets there are, this gives how many bodies there are in the system
    
        
    #loop through all points
    for k in range(0,count):
        st=f.index("{",start_count_index) #Get first body information
        en=f.index("}",end_count_index) #Stop body information once "}" is found
    
        #Create list of each body's information
        x=[]
        for i in range(st+1,en):
            y=f[i].partition(":")
            x.append({y[0].strip("\'").strip().strip(" ' "):y[2].strip()})        

    for i in x:
        if  list(i.keys())[0] == 'ft':
            ft=(list(i.values())[0])
        if  list(i.keys())[0] == 'df':
            dft=(list(i.values())[0])        
        if  list(i.keys())[0] == 'ddf':
            ddft=(list(i.values())[0])                   
        if list(i.keys())[0] == 'function':
            print("read function")
            ft_file=os.getcwd()+"\\f_t.txt"
            
            f=open(ft_file,'r')
            contents=f.read()
            f.close()                     
                
           
    return float(ft),float(dft),float(ddft)