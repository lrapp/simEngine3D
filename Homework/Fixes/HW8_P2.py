from simEngine3D import *
import time as ttime
import matplotlib.pyplot as plt
import numpy as np
import os 



t_start=0
t_end=2
t_step=0.001


file_dir=os.getcwd()
file_name="\\Homework\\Fixes\\HW8_P2_input.txt"
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

# output_file=os.getcwd()+"\\out.csv"
# SYS.outputs.to_csv(output_file)
#%%

position=SYS.outputs["r"]
x=[]
y=[]
z=[]
x1=[]
y1=[]
z1=[]
for i in position:
    x.append(i[0])
    y.append(i[1])
    z.append(i[2])
    x1.append(i[3])
    y1.append(i[4])
    z1.append(i[5])    
    
    
    
time_list=np.arange(t_start,t_end,t_step)
plt.plot(time_list,z1)



t=time_list[0::50]
yyy=y[0::50]
zzz=z[0::50]

y2=y1[0::50]
z2=z1[0::50]

plt.figure()
dist=[]
for i in range(0,len(t)):
    plt.xlim(-6,6)
    plt.ylim(-6,6)
    yy=[0,yyy[i]*2]
    zz=[0,zzz[i]*2]
    
    yy2=[yyy[i]*2,y2[i]]
    zz2=[zzz[i]*2,z2[i]]
    
    line_1=plt.plot(yy,zz,c='b')
    line_2=plt.plot(yy2,zz2,c='r')
    # dist.append(np.sqrt((yy[1]-yy2[1])**2+(zz[1]-zz2[1])**2))
    scat_1=plt.plot(yyy[i],zzz[i],'o',c='b')
    scat_2=plt.plot(y2[i],z2[i],'o',c='r')    
    plt.pause(0.1)
    line=line_1.pop()
    line.remove()
    line=line_2.pop()
    line.remove()
    #%%
#plot x,y,z for body 2 - position vs time
#plt.plot(time_list,x1)
# plt.plot(time_list,y1)
# plt.plot(time_list,z1)
    
