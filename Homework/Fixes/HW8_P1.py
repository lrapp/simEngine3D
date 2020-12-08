import os 

os.chdir("..")
os.chdir("..")

from simEngine3D import *
import time as ttime
import matplotlib.pyplot as plt
import numpy as np




t_start=0 
t_end=10
t_step=0.01


file_dir=os.getcwd()
file_name="\\Homework\\Fixes\\HW8_P1_input.txt"
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



position=SYS.outputs["r"]
x=[]
y=[]
z=[]
for i in position:
    x.append(i[0])
    y.append(i[1])
    z.append(i[2])
    
    
    
time_list=np.arange(t_start,t_end,t_step)
plt.plot(time_list,y)