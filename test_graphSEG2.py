# -*- coding: utf-8 -*-
"""
Created on Sun Dec  6 16:28:49 2020

@author: tonianr
"""
#%% Directory
import os 

dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)



#%% Plot
dat_file = "./exampledata/302.dat"
from graphSEG2 import SEG2grapher
SEG2grapher.plot(dat_file,15)
#%% Pick
dat_file = "./exampledata_FGSeismic2019/108.dat"
from graphSEG2 import SEG2grapher
SEG2grapher.pick(dat_file,25,'false')
#%% Analyze
from graphSEG2 import SEG2grapher
SEG2grapher.analyze(pckfile = "exampledata/303.txt",show_autocrossover='true')
#%% Crossover
pck_file = "exampledata_FGSeismic2019/109.txt"
from graphSEG2 import SEG2grapher
SEG2grapher.crossover(pck_file, show_autocrossover='true')
#%% Closestpoint

from findclosestpoint import closestpoint
import numpy as np

x = [1,2,3,4,5,6,8]
y = [1,2,3,4,5,6,8]
a = np.stack((x,y))
a = a.T
pt = [3.12, 3.12]

w = closestpoint(pt,a)
print(w)


