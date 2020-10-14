

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
#    if e0 != 0:       
#        e1=(A[2,1]-A[1,2])/(4*e0)
#        e2=(A[0,2]-A[2,0])/(4*e0)
#        e3=(A[1,0]-A[0,1])/(4*e0)
#    elif (1+2*A[0,0]-np.trace(A))/4 != 0:
#        e1=np.sqrt((1+2*A[0,0]-np.trace(A))/4)
#        e2=(A[1,0]+A[0,1])/(4*e1)
#        e3=(A[2,0]+A[0,2])/(4*e1)
#    elif (1+2*A[1,1]-np.trace(A))/4 !=0:
#        e2=np.sqrt((1+2*A[1,1]-np.trace(A))/4)
#        e1=(A[1,0]+A[0,1])/(4*e2)
#        e3=(A[2,1]+A[1,2])/(4*e2)
#    elif (1+2*A[2,2]-np.trace(A))/4 != 0:
#        e3=np.sqrt((1+2*A[2,2]-np.trace(A))/4)
#        e1=(A[2,0]+A[0,2])/(4*e3)
#        e2=(A[2,1]+A[1,2])/(4*e3)
             
    p=np.empty([4,1])
    p[:,0]=np.array([e0,e1,e2,e3])

    return p

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

