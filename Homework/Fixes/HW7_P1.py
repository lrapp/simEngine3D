import sys
sys.path.append('../..')
import time as ttime
import matplotlib.pyplot as plt
import numpy as np
import os
from simEngine3D import *


if any('SPYDER' in name for name in os.environ):
    print("in Spyder, inputs are entered here, not in command line")
    t_start=0
    t_end=1
    t_step=0.001    
    
else:
    print("enter end time of simulation and hit enter. \nExample: 10  (takes ~266 seconds)")
    t_end=float(input())
    print("enter stepsize and hit enter. \nExample: 0.001")
    t_step=float(input())

print("\n")
print("inputs:","t_end=",str(t_end),",","t_step=",str(t_step))


t_start=0

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
plt.show()
#%%
position=SYS.outputs.r
x=[]
y=[]
z=[]
for i in range(0,len(SYS.outputs)):
    x.append(position[i][0])    
    y.append(position[i][1])
    z.append(position[i][2])    
    
    
#this checks to see if file is running in Spyder or from command line
#if in command line, it will wait for input from user before closing figures
if any('SPYDER' in name for name in os.environ):
    print("in Spyder, don't need to hold figures open")
else:
    print("press enter to continue")
    input()

           