import sys
sys.path.append('../..')
import os
from simEngine3D import *
import time as ttime
import matplotlib.pyplot as plt
import numpy as np






if any('SPYDER' in name for name in os.environ):
    print("in Spyder, inputs are entered here, not in command line")
    t_start=0
    t_end=10
    t_step=0.01   
    
else:
    print("enter end time of simulation and hit enter. \nExample: 10  (takes ~16.5 seconds)")
    t_end=float(input())
    print("enter stepsize and hit enter. \nExample: 0.01")
    t_step=float(input())

print("\n")
print("inputs:","t_end=",str(t_end),",","t_step=",str(t_step))
t_start=0      


file_dir=os.getcwd()
file_name="\\HW8_P1_input.txt"
file=file_dir+file_name

SYS=sys() #Create system object
#%% set derivites of driving function specific to this problem
def df(self):
    return 0 #no driving function here
    

def ddf(self):
    return 0 #no driving function here
 
setattr(SYS,'df',df)     
setattr(SYS,'ddf',ddf)         
#%%
tic=ttime.perf_counter()
SYS=dynamic_analysis(SYS,t_start,t_step,t_end,file)
toc=ttime.perf_counter()

elapsed_time=toc-tic
print("time elapsed=",elapsed_time)

#%% Plot positions
time_list=np.arange(t_start,t_end,t_step)

position=SYS.outputs["r"]
x=[]
y=[]
z=[]
for i in position:
    x.append(i[0])
    y.append(i[1])
    z.append(i[2])
    
plt.figure()    
plt.xlim(0,10)
plt.title("HW8_P1 O' Position")
plt.xlabel("Time (s)")
plt.ylabel("Position")
plt.plot(time_list,x,label="x")
plt.plot(time_list,y,label="y")
plt.plot(time_list,z,label="z")
plt.legend()
    

#%%    Plot velocities
dr=SYS.outputs["dr"]

dx=[]
dy=[]
dz=[]
for i in dr:
    dx.append(i[0])
    dy.append(i[1])
    dz.append(i[2])
    
plt.figure()    
plt.xlim(0,10)
plt.title("HW8_P1 O' velocities")
plt.xlabel("Time (s)")
plt.ylabel("Omega (rad/s)")
plt.plot(time_list,dx,label="dx")
plt.plot(time_list,dy,label="dy")
plt.plot(time_list,dz,label="dz")
plt.legend()
    
#%%
ddr=SYS.outputs["ddr"]
ddx=[]
ddy=[]
ddz=[]
for i in ddr:
    ddx.append(i[0])
    ddy.append(i[1])
    ddz.append(i[2])

plt.figure()    
plt.xlim(0,10)
plt.title("HW8_P1 O' Accelerations")
plt.xlabel("Time (s)")
plt.ylabel("Acceleration (rad/s^2)")
plt.plot(time_list,ddx,label="ddx")
plt.plot(time_list,ddy,label="ddy")
plt.plot(time_list,ddz,label="ddz")
plt.legend()


#%%


if t_step <= 0.01:
    spare_index=10
else:
    spare_index=5
    
t=time_list[0::spare_index]
y_1_sparse=y[0::spare_index]
z_1_sparse=z[0::spare_index]


plt.figure()
for i in range(0,len(t)):
    plt.xlim(-6,6)
    plt.ylim(-6,6)
    yy=[0,y_1_sparse[i]*2]
    zz=[0,z_1_sparse[i]*2]
    
    line_1=plt.plot(yy,zz,c='b')
    scat_1=plt.plot(y_1_sparse[i],z_1_sparse[i],'o',c='b')

    plt.pause(0.01)
    plt.title("time="+str(round(t[i],1)))
    line=line_1.pop()
    line.remove()    
    
    
#this checks to see if file is running in Spyder or from command line
#if in command line, it will wait for input from user before closing figures
if any('SPYDER' in name for name in os.environ):
    print("in Spyder, don't need to hold figures open")
else:
    print("press enter to continue")
    input()

