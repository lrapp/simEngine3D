

def build_p(theta):
    import numpy as np
    A=np.zeros([3,3])
    A[0,2]=-1
    A[1,0]=np.sin(theta)
    A[1,1]=-np.sin(theta)
    A[2,0]=-np.cos(theta)
    A[2,1]=-np.sin(theta)
    
    e0=np.sqrt((np.trace(A)+1)/4)
    e1=(A[2,1]-A[1,2])/(4*e0)
    e2=(A[0,2]-A[2,0])/(4*e0)
    e3=(A[1,0]-A[0,1])/(4*e0)
    p=np.empty([4,1])
    p[:,0]=np.array([e0,e1,e2,e3])
    return p

