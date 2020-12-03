from simEngine3D import *
import time as ttime
import matplotlib.pyplot as plt
import numpy as np




t_start=0
t_end=10
t_step=0.1

tic=ttime.perf_counter()
SYS=dynamic_analysis(t_start,t_step,t_end)
toc=ttime.perf_counter()

elapsed_time=toc-tic
print("time elapsed=",elapsed_time)



position=SYS.outputs["r0"]
x=[]
y=[]
z=[]
for i in position:
    x.append(i[0])
    y.append(i[1])
    z.append(i[2])
    
    
    
time_list=np.arange(t_start,t_end,t_step)
plt.plot(time_list,y)
