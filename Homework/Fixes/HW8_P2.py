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
file_name="\\HW8_P2_input.txt"
file=file_dir+file_name

# file="C:\\Users\\Logan\\Desktop\\simEngine3D\\HW8_P1_revJoint.txt"
SYS=sys()
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

#%%
#Output can be saved easily to .csv, pickle, or your favorite file type.
# output_csv_file=os.getcwd()+"\\out.csv"
# output_pickle_file=os.getcwd()+"\\out.pkl"
# SYS.outputs.to_csv(output_csv_file)
# SYS.outputs.to_pickle(output_pickle_file)
#%%

#Get positions from output DataFrame
tic=ttime.perf_counter()
position=SYS.outputs["r"].to_numpy()
p_a=np.concatenate(position)

x1=p_a[0::6]
y1=p_a[1::6]
z1=p_a[2::6]

x2=p_a[3::6]
y2=p_a[4::6]
z2=p_a[5::6]

toc=ttime.perf_counter()

# print("array time=",toc-tic)

#The above "array" method is faster than the loop method below, so I've commented it out.
# tic=ttime.perf_counter()
# x=[]
# y=[]
# z=[]
# x1=[]
# y1=[]
# z1=[]
# for i in position:
#     x.append(i[0])
#     y.append(i[1])
#     z.append(i[2])
#     x1.append(i[3])
#     y1.append(i[4])
#     z1.append(i[5])    
    
# toc=ttime.perf_counter()


# print("list time=",toc-tic)   

#%%
time_list=np.arange(t_start,t_end,t_step)


position=SYS.outputs["r"]
x=[]
y=[]
z=[]

for i in position:
    x.append(i[0::3])
    y.append(i[1::3])
    z.append(i[2::3])
    

x_array=np.concatenate(x)
y_array=np.concatenate(y)
z_array=np.concatenate(z)
x_list=[]
y_list=[]
z_list=[]
for i in range(0,SYS.nb):
    x_list.append(x_array[i::SYS.nb])
    y_list.append(y_array[i::SYS.nb])
    z_list.append(z_array[i::SYS.nb])



    
fig=plt.figure(figsize=(15,10))   
ax1=fig.add_subplot(3,2,1) 
ax1.set(xlim=(0,10))
ax1.set_title("HW8_P2 O1' Position")
ax1.set(ylabel="Position")
ax1.plot(time_list,x_list[0],label="x1")
ax1.plot(time_list,y_list[0],label="y1")
ax1.plot(time_list,z_list[0],label="z1")
ax1.legend()
    
ax2=fig.add_subplot(3,2,2)
ax2.set(xlim=(0,10))
ax2.set_title("HW8_P2 O2' Position")
ax2.set(ylabel="Position")
ax2.plot(time_list,x_list[1],label="x2")
ax2.plot(time_list,y_list[1],label="y2")
ax2.plot(time_list,z_list[1],label="z2")
ax2.legend()
    

#    Plot velocities
dr=SYS.outputs["dr"]

dx=[]
dy=[]
dz=[]
for i in dr:
    dx.append(i[0::3])
    dy.append(i[1::3])
    dz.append(i[2::3])


dx_array=np.concatenate(dx)
dy_array=np.concatenate(dy)
dz_array=np.concatenate(dz)
dx_list=[]
dy_list=[]
dz_list=[]
for i in range(0,SYS.nb):
    dx_list.append(dx_array[i::SYS.nb])
    dy_list.append(dy_array[i::SYS.nb])
    dz_list.append(dz_array[i::SYS.nb])
    
    

ax1=fig.add_subplot(3,2,3) 
ax1.set(xlim=(0,10))
ax1.set_title("HW8_P1 O1' velocities")
ax1.set(ylabel="Omega (rad/s)")
ax1.plot(time_list,dx_list[0],label="dx1")
ax1.plot(time_list,dy_list[0],label="dy1")
ax1.plot(time_list,dz_list[0],label="dz1")
ax1.legend()
    
ax2=fig.add_subplot(3,2,4)
ax2.set(xlim=(0,10))
ax2.set_title("HW8_P1 O2' velocities")
ax2.set(ylabel="Omega (rad/s)")
ax2.plot(time_list,dx_list[1],label="dx2")
ax2.plot(time_list,dy_list[1],label="dy2")
ax2.plot(time_list,dz_list[1],label="dz2")
ax2.legend()
    
    
# plt.figure()    
# plt.xlim(0,10)
# plt.title("HW8_P1 O1' velocities")
# plt.xlabel("Time (s)")
# plt.ylabel("Omega (rad/s)")
# index=0
# plt.plot(time_list,dx_list[index],label="dx")
# plt.plot(time_list,dy_list[index],label="dy")
# plt.plot(time_list,dz_list[index],label="dz")
# plt.legend()
    

# plt.figure()    
# plt.xlim(0,10)
# plt.title("HW8_P1 O2' velocities")
# plt.xlabel("Time (s)")
# plt.ylabel("Omega (rad/s)")
# index=1
# plt.plot(time_list,dx_list[index],label="dx")
# plt.plot(time_list,dy_list[index],label="dy")
# plt.plot(time_list,dz_list[index],label="dz")
# plt.legend()
    

ddr=SYS.outputs["ddr"]

ddx=[]
ddy=[]
ddz=[]
for i in ddr:
    ddx.append(i[0::3])
    ddy.append(i[1::3])
    ddz.append(i[2::3])


ddx_array=np.concatenate(ddx)
ddy_array=np.concatenate(ddy)
ddz_array=np.concatenate(ddz)
ddx_list=[]
ddy_list=[]
ddz_list=[]
for i in range(0,SYS.nb):
    ddx_list.append(ddx_array[i::SYS.nb])
    ddy_list.append(ddy_array[i::SYS.nb])
    ddz_list.append(ddz_array[i::SYS.nb])
    
    
ax1=fig.add_subplot(3,2,5) 
ax1.set(xlim=(0,10))
ax1.set_title("HW8_P1 O1' Acceleration")
ax1.set(xlabel="Time (s)",ylabel="Accleration (rad/s^2)")
ax1.plot(time_list,ddx_list[0],label="ddx1")
ax1.plot(time_list,ddy_list[0],label="ddy1")
ax1.plot(time_list,ddz_list[0],label="ddz1")
ax1.legend()
    
ax2=fig.add_subplot(3,2,6)
ax2.set(xlim=(0,10))
ax2.set_title("HW8_P1 O2' Acceleration")
ax2.set(xlabel="Time (s)", ylabel="Accleration (rad/s^2)")
ax2.plot(time_list,ddx_list[1],label="ddx2")
ax2.plot(time_list,ddy_list[1],label="ddy2")
ax2.plot(time_list,ddz_list[1],label="ddz2")
ax2.legend()
        

# plt.figure()    
# plt.xlim(0,10)
# plt.title("HW8_P1 O1' Acceleration")
# plt.xlabel("Time (s)")
# plt.ylabel("Accleration (rad/s^2)")
# index=0
# plt.plot(time_list,ddx_list[index],label="ddx")
# plt.plot(time_list,ddy_list[index],label="ddy")
# plt.plot(time_list,ddz_list[index],label="ddz")
# plt.legend()
    

# plt.figure()    
# plt.xlim(0,10)
# plt.title("HW8_P1 O2' Acceleration")
# plt.xlabel("Time (s)")
# plt.ylabel("Accleration (rad/s^2)")
# index=1
# plt.plot(time_list,ddx_list[index],label="ddx")
# plt.plot(time_list,ddy_list[index],label="ddy")
# plt.plot(time_list,ddz_list[index],label="ddz")
# plt.legend()
    
    #%%
time_list=np.arange(t_start,t_end,t_step)
# plt.plot(time_list,z1)

y1=y_list[0]
z1=z_list[0]

y2=y_list[1]
z2=z_list[1]

spare_index=50
t=time_list[0::spare_index]
y_1_sparse=y1[0::spare_index]
z_1_sparse=z1[0::spare_index]

y_2_sparse=y2[0::spare_index]
z_2_sparse=z2[0::spare_index]

fig=plt.figure()
plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
plt.xlabel('y position',size=20)
plt.ylabel('z position',size=20)
dist=[]

for i in range(0,len(t)):
    plt.xlim(-6,6)
    plt.ylim(-6,6)
    plt.title("time="+str(round(t[i],1)))
    yy=[0,y_1_sparse[i]*2]
    zz=[0,z_1_sparse[i]*2]
    
    # yy2=[y_1_sparse[i]*2,y_2_sparse[i]]
    # zz2=[z_1_sparse[i]*2,z_2_sparse[i]]
    
    yy2=[y_1_sparse[i]*2,y_2_sparse[i]]
    zz2=[z_1_sparse[i]*2,z_2_sparse[i]]
    
    line_1=plt.plot(yy,zz,c='b')
    line_2=plt.plot(yy2,zz2,c='r')
    # dist.append(np.sqrt((yy[1]-yy2[1])**2+(zz[1]-zz2[1])**2))
    scat_1=plt.plot(y_1_sparse[i],z_1_sparse[i],'o',c='b')
    scat_2=plt.plot(y_2_sparse[i],z_2_sparse[i],'o',c='r')    
    plt.pause(0.01)
    
    if i < len(t)-1:
        line=line_1.pop()
        line.remove()
        line=line_2.pop()
        line.remove()

#this checks to see if file is running in Spyder or from command line
#if in command line, it will wait for input from user before closing figures
if any('SPYDER' in name for name in os.environ):
    print("in Spyder, don't need to hold figures open")
else:
    print("press enter to continue")
    input()

       