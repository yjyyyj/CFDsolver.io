#!/bin/usr/python3
import numpy as np
import matplotlib.pyplot as plt
import pathlib
import pprint
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import Normalize
import ffmpeg as fp

###### plot set ##################################################

range_min=-1e-4
range_max= 1e-4

fname = "output_000000.dat"
files = np.sort(list(pathlib.Path('./').glob('output_*.dat')))
fnames = [ i.name for i in files]
print(fnames)
# global settings
plt.rcParams['xtick.direction'] = 'in'#x軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams['ytick.direction'] = 'in'#y軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams['xtick.major.width'] = 1.0 #x軸主目盛り線の線幅
plt.rcParams['ytick.major.width'] = 1.0 #y軸主目盛り線の線幅
plt.rcParams['font.size'] = 22 #フォントの大きさ
plt.rcParams['axes.linewidth'] = 1.0 # 軸の線幅edge linewidth。囲みの太さ
plt.rcParams['text.usetex'] = True # Latex text
plt.rcParams["axes.formatter.use_mathtext"]=True

# *** read init data ***********************
data0=(np.loadtxt(fname,dtype="float",skiprows=0)).T
print(np.shape(data0))

# data definition
jmax = 501
kmax = 3
resid = np.zeros((len(files),len(data0)+1))

# *** loop start *************************
for i, fname in enumerate(fnames):

  # data read
  data1=(np.loadtxt(fname,dtype="float",skiprows=0)).T

  x1=data1[0,:jmax] # x
  x2=data1[0,:kmax]  # y
  diff = data1 -data0  # diff
  var = np.average(diff*diff,axis=1)
  
  resid[i] = np.append(i,np.sqrt(var))

# *****************************************************

# write norm file
np.savetxt('norm.dat', resid)
