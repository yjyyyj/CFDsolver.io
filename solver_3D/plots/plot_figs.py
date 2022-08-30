#!/bin/usr/python3
import numpy as np
import matplotlib.pyplot as plt
import pathlib
import pprint
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import Normalize
import ffmpeg as fp

###### make GIF ##################################################

stream = fp.input("./img%03d.png", framerate=10)
stream = fp.output(stream, "anime.mp4", pix_fmt='yuv420p')

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


# *** loop start *************************
for i, fname in enumerate(fnames):
  fig = plt.figure()
  fig.set_size_inches([12, 10])

  plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.16)

  # data read
  data1=(np.loadtxt(fname,dtype="float",skiprows=0)).T
  print(np.shape(data1))

  # data definition
  jmax = 401
  kmax = 401

  x1=data1[0,:jmax]  # x
  x2=x1             # y
  q1=data1[2]  # r
  q2=data1[3]  # u
  q3=data1[4]  # v
  q4=data1[5]-1  # p
  q5=data1[6]  # ry1
  q6=data1[7]  # ry2
  q7=data1[8]  # ry3

  # *** data read for diff *****************************
  # data3=(np.loadtxt(fname3,dtype="float",skiprows=0)).T
  # data4=(np.loadtxt(fname4,dtype="float",skiprows=0)).T
  # data5=(np.loadtxt(fname5,dtype="float",skiprows=0)).T

  # x1=data3[0]  # x
  # # p_diff
  # q1=data3[3]  # div
  # q2=data4[3]  # p
  # q3=data5[3]  # proposed

  # *****************************************************

  q1 = np.array(q1).reshape(-1, jmax).tolist()
  q2 = np.array(q2).reshape(-1, jmax).tolist()
  q3 = np.array(q3).reshape(-1, jmax).tolist()
  q4 = np.array(q4).reshape(-1, jmax).tolist()
  q5 = np.array(q5).reshape(-1, jmax).tolist()
  q6 = np.array(q6).reshape(-1, jmax).tolist()
  q7 = np.array(q7).reshape(-1, jmax).tolist()


  # format
  plt.xlim(0,1)
  # plt.ylim(range_min,range_max)
  plt.grid(False)
  # nmax = x1.shape[0]-6
  nmax = x1.shape[0]

  plt.title(r'$t=$'+ "{:.2f}".format(0.1*i))
  # plt.title("Plot 2D array")
  plt.xlabel(r'$x$',fontsize=28,y=0)
  plt.ylabel(r"$y$",fontsize=28)
  # plt.xticks(np.arange(0,1.2,0.2))

  xx, yy = np.meshgrid(x1, x2)
  cf = plt.contourf(xx, yy, q1, cmap="jet")
  # cf = plt.contourf(xx, yy, q2, cmap="jet", levels=12)
  # cf = plt.contourf(xx, yy, q3, cmap="jet", levels=12)
  # cf = plt.contourf(xx, yy, q4, cmap="jet", levels=12)
  # cf = plt.contourf(xx, yy, q5, cmap="jet", levels=12)
  # cf = plt.contourf(xx, yy, q6, cmap="jet", levels=12)
  # cf = plt.contourf(xx, yy, q7, cmap="jet", levels=12)

  ypoint = int((kmax-1)/2)
  plt.plot(x1[:nmax], q1[ypoint], linewidth=2, color="blue", linestyle="solid",zorder=1)
  plt.plot(x1[:nmax], q5[ypoint], linewidth=2, color="green", linestyle="solid",zorder=1)
  plt.plot(x1[:nmax], q6[ypoint], linewidth=2, color="pink", linestyle="solid",zorder=1)
  plt.plot(x1[:nmax], q7[ypoint], linewidth=2, color="orange", linestyle="solid",zorder=1)

  cb = plt.colorbar(cf)
  cb.set_label(r'$\rho$',fontsize=28,y=1,x=1.2,rotation=0)
  # cb.set_label(r"$(p-p_{\mathrm{exact}})/p_{\mathrm{exact}}$",fontsize=28)
  # cb.set_label(r"$(u-u_{\mathrm{exact}})/u_{\mathrm{exact}}$",fontsize=28)
  cb.formatter.set_scientific(True)


  # plt.show()
  fig.savefig("img{:03d}.png".format(i))

fp.run(stream)

