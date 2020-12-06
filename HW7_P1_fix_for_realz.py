from simEngine3D import *
import time as ttime
import matplotlib.pyplot as plt
import numpy as np
import os 



t_start=0
t_end=5
t_step=0.001

file_dir=os.getcwd()
file_name="\\Homework\\Fixes\\revJoint_fix2.txt"
file=file_dir+file_name

# file="C:\\Users\\Logan\\Desktop\\simEngine3D\\Homework\\Fixes\\revJoint_fix2.txt"

#create system object
SYS=sys()
#%%
def df(t):
    return -np.pi/2. * np.sin(2.* t) * np.cos(np.pi/4.* np.cos(2.* t))

    

def ddf(t):
    return -np.pi/4. * (4.* np.cos(2.* t) * np.cos(np.pi/4.* np.cos(2* t)) + \
            np.pi * np.sin(2.* t) * np.sin(2.* t) * np.sin(np.pi/4.* np.cos(2.* t)))

 
setattr(SYS,'df',df)     
setattr(SYS,'ddf',ddf)         
#%%
tic=ttime.perf_counter()
SYS=inverse_dynamics_analysis(SYS,t_start,t_step,t_end,file)
toc=ttime.perf_counter()

elapsed_time=toc-tic
print("time elapsed=",elapsed_time)


torque=SYS.torque

time_list=np.arange(t_start,t_end,t_step)
plt.plot(time_list,torque[:,0])

#%%
y=[]
z=[]
for i in range(0,len(out)):
    y.append(out.r0[i][1][0])
    z.append(out.r0[i][2][0])    