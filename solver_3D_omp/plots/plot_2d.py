#!/bin/usr/python3
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import Normalize # Normalizeをimport

range_min=-1e-2
range_max= 1e-3
fname1 = "output_000000.dat"
# fname1 = "output_010000.dat"
# fname1 = "output_pro.dat"
# fname1 = "output_conv.dat"


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
fig.set_size_inches([10, 8])
# fig.figsize= fig.figaspect(2)
plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.16)

# data read
data1=(np.loadtxt(fname1,dtype="float",skiprows=0)).T
print(np.shape(data1))

# data definition
jmax = 101
kmax = 101
# jmax = 501
# kmax = 501

x1=data1[0,:jmax]  # x
x2=data1[0,:kmax]  # y
q1=data1[2]  # r
q2=data1[3]  # u
q3=data1[4]  # v
q4=data1[5]  # w
q5=data1[6]  -1# p
q6=data1[7]  # ry1
q7=data1[7]  # ry2
# q7=data1[8]  # ry3

q1 = np.array(q1).reshape(-1, jmax).tolist()
q2 = np.array(q2).reshape(-1, jmax).tolist()
q3 = np.array(q3).reshape(-1, jmax).tolist()
q4 = np.array(q4).reshape(-1, jmax).tolist()
q5 = np.array(q5).reshape(-1, jmax).tolist()
q6 = np.array(q6).reshape(-1, jmax).tolist()
q7 = np.array(q7).reshape(-1, jmax).tolist()

# format
# plt.xlim(0,1)
plt.grid(False)
nmax = x1.shape[0]

# plt.title("Plot 2D array")
plt.xlabel(r'$x$',fontsize=22)
plt.ylabel(r"$y$",fontsize=22)
# plt.xticks(np.arange(0,1.2,0.2))
plt.xticks(y=-0.02)

xx, yy = np.meshgrid(x1, x2)
levels = [0.06,0.1199,0.18,0.24,0.30,0.36,0.42,0.48,0.54,0.60,0.66,0.72,0.78,0.8401,0.90]
cf = plt.contourf(xx, yy, q1, levels, cmap="jet")

# cf = plt.contourf(xx, yy, q1, cmap="jet", levels=14)
# cf = plt.contourf(xx, yy, q2, cmap="jet", levels=12)
# cf = plt.contourf(xx, yy, q3, cmap="jet", levels=12)
# cf = plt.contourf(xx, yy, q4, cmap="jet", levels=12)
# cf = plt.contourf(xx, yy, q5, cmap="jet", levels=12)
# cf = plt.contourf(xx, yy, q6, cmap="jet", levels=12)

# cf = plt.pcolormesh(xx, yy, q1, cmap="jet", shading="gouraud")
# cf = plt.pcolormesh(xx, yy, q2, cmap="bwr", shading="gouraud")
# cf = plt.pcolormesh(xx, yy, q3, cmap="bwr", shading="gouraud")
# cf = plt.pcolormesh(xx, yy, q4, cmap="bwr", shading="gouraud")
# cf = plt.pcolormesh(xx, yy, q5, cmap="bwr", shading="gouraud")
# cf = plt.pcolormesh(xx, yy, q6, cmap="bwr", shading="gouraud")

cb = fig.colorbar(cf)
cb.set_label(r'$\rho$',fontsize=22,loc='center', rotation=270,labelpad=40)
# cb.set_label(r'$(p-p_{\rm{exact}})/p_{\rm{exact}}$',fontsize=24,loc='center', rotation=270,labelpad=20)

cb.formatter.set_scientific(True)
# cb.ax.ticklabel_format(style='sci', scilimits=(-6,6))
# cb.set_clim(0.11, 0.85) # カラースケールのグラデーションの端点
# plt.clim(0.06, 0.9) # カラースケールのグラデーションの端点
# plt.clim(-5e-12, 5e-12) # カラースケールのグラデーションの端点


ax = plt.gca()
aspect = 1/1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])                     
ax.set_aspect(aspect)

# plt.show()
fig.savefig("2dtest.png")
