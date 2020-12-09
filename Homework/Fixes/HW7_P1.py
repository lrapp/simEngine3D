import sys
sys.path.append('../..')
import time as ttime
import matplotlib.pyplot as plt
import numpy as np
import os
from simEngine3D import *




t_start=0
t_end=10
t_step=0.001

file_dir=os.getcwd()
file_name="\\revJoint_fix2.txt"
file=file_dir+file_name


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
print("Solving...")
SYS=inverse_dynamics_analysis(SYS,t_start,t_step,t_end,file)
toc=ttime.perf_counter()

elapsed_time=toc-tic
print("time elapsed=",elapsed_time)


torque=SYS.torque

time_list=np.arange(t_start,t_end,t_step)
fig=plt.figure()
plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
plt.xlabel('Time (s)',size=20)
plt.ylabel('Reaction Torque (Nm)',size=20)
plt.plot(time_list,torque[:,0],label='x')
plt.plot(time_list,torque[:,1],label='y')
plt.plot(time_list,torque[:,2],label='z')
plt.legend()
#%%
position=SYS.outputs.r
x=[]
y=[]
z=[]
for i in range(0,len(SYS.outputs)):
    x.append(position[i][0])    
    y.append(position[i][1])
    z.append(position[i][2])    