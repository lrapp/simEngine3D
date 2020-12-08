import os
os.chdir("..")
os.chdir("..")
from simEngine3D import *
import time as ttime
import matplotlib.pyplot as plt
import numpy as np



t_start=0
t_end=5
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
# plt.plot(time_list,z1)



t=time_list[0::50]
y_1_sparse=y1[0::50]
z_1_sparse=z1[0::50]

y_2_sparse=y2[0::50]
z_2_sparse=z2[0::50]

plt.figure()
dist=[]
for i in range(0,len(t)):
    plt.xlim(-6,6)
    plt.ylim(-6,6)
    yy=[0,y_1_sparse[i]*2]
    zz=[0,z_1_sparse[i]*2]
    
    yy2=[y_1_sparse[i]*2,y_2_sparse[i]]
    zz2=[z_1_sparse[i]*2,z_2_sparse[i]]
    
    line_1=plt.plot(yy,zz,c='b')
    line_2=plt.plot(yy2,zz2,c='r')
    # dist.append(np.sqrt((yy[1]-yy2[1])**2+(zz[1]-zz2[1])**2))
    scat_1=plt.plot(y_1_sparse[i],z_1_sparse[i],'o',c='b')
    scat_2=plt.plot(y_2_sparse[i],z_2_sparse[i],'o',c='r')    
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
    
