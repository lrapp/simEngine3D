# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 14:27:42 2020
Restarting with pieces of my SimEngine3D from homework.  Starting over here for a clean slate.

@author: Logan
"""
import os
os.chdir("C:\\Users\\Logan\\Desktop\\simEngine3D")
from simEngine3D_dataload import data_file, constraints_in
import tkinter as tk
from tkinter import ttk



X=data_file("C:\\Users\\Logan\\Desktop\\simEngine3D\\revJoint_fix.txt") #list of body objects
constraint_list=constraints_in("C:\\Users\\Logan\\Desktop\\simEngine3D\\revJoint_fix.txt") #list of constraints

