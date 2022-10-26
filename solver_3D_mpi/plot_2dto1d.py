#!/bin/usr/python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

range_min=-1e-2
range_max= 1e-3
fname1 = "output_000000.dat"
# fname1 = "output_000005.dat"
# fname1 = "output_006000.dat"
# fname1 = "output_008000.dat"
# fname1 = "output_050000.dat"
# fname1 = "output_010000.dat"


# global settings
plt.rcParams['xtick.direction'] = 'in'#x軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams['ytick.direction'] = 'in'#y軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams['xtick.major.width'] = 1.0 #x軸主目盛り線の線幅
plt.rcParams['ytick.major.width'] = 1.0 #y軸主目盛り線の線幅
plt.rcParams['font.size'] = 22 #フォントの大きさ
plt.rcParams['axes.linewidth'] = 1.0 # 軸の線幅edge linewidth。囲みの太さ
plt.rcParams['text.usetex'] = True # Latex text
plt.rcParams["axes.formatter.use_mathtext"]=True

fig = plt.figure()
fig.set_size_inches([12, 10])
plt.subplots_adjust(left=0.18, right=0.95, top=0.9, bottom=0.16)

# data read
# data0=(np.loadtxt(fname0,dtype="float",skiprows=0)).T
data1=(np.loadtxt(fname1,dtype="float",skiprows=0)).T
print(np.shape(data1))

# data definition
jmax = 501
kmax = 3

x1=data1[0,:jmax]  # x

q1=data1[2]  # r
q2=data1[3]  # u
q3=data1[4]  # v
q4=data1[5]  # w
q5=data1[6]  # p
q6=data1[7]  # ry1
q7=data1[8]  # ry2

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

# plt.title("Plot 2D array")
plt.xlabel(r'$x$',fontsize=22,y=0)
# plt.ylabel(r"$y$",fontsize=22)
# plt.xticks(np.arange(0,1.2,0.2))

ypoint = int((kmax-1)/2)
plt.plot(x1[:nmax], q1[ypoint], linewidth=2, color="blue", linestyle="solid",zorder=1)
plt.plot(x1[:nmax], q2[ypoint], linewidth=2, color="red", linestyle="solid",zorder=1)
plt.plot(x1[:nmax], q4[ypoint], linewidth=2, color="orange", linestyle="solid",zorder=1)
plt.plot(x1[:nmax], q5[ypoint], linewidth=2, color="green", linestyle="solid",zorder=1)
plt.plot(x1[:nmax], q6[ypoint], linewidth=2, color="cyan", linestyle="solid",zorder=1)
plt.plot(x1[:nmax], q7[ypoint], linewidth=2, color="pink", linestyle="solid",zorder=1)

# plt.plot(x1[:nmax], data0[6,:jmax], "o", color="black", markerfacecolor="white",zorder=-1)
# plt.plot(x1[:nmax], data0[7,:jmax], "s", color="black", markerfacecolor="white",zorder=-1)

####### 2D contour ###############################################
# xx, yy = np.meshgrid(x1, x2)
# cf = plt.contourf(xx, yy, q1, cmap="jet", levels=18)
# # cf = plt.contourf(xx, yy, q2, cmap="jet", levels=12)
# # cf = plt.contourf(xx, yy, q3, cmap="jet", levels=12)
# # cf = plt.contourf(xx, yy, q4, cmap="jet", levels=12)
# # cf = plt.contourf(xx, yy, q5, cmap="jet", levels=12)
# # cf = plt.contourf(xx, yy, q6, cmap="jet", levels=12)
# # cf = plt.contourf(xx, yy, q7, cmap="jet", levels=12)


# cb = plt.colorbar(cf)
# cb.set_label(r'$\rho$',fontsize=22,y=1,x=1.2,rotation=0)
# cb.formatter.set_scientific(True)
# # cb.ax.ticklabel_format(style='sci', scilimits=(-6,6)) 
# # plt.show()

fig.savefig("2dtest_1d.png")
