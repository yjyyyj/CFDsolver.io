#!/bin/usr/python3
from json.encoder import INFINITY
from re import A
from weakref import ref
import matplotlib
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


fname1 = "output_000000.dat"
# fname1 = "output_001000.dat"
# fname1 = "output_020000.dat"
fname1 = "output_060000.dat"
# fname1 = "output_200000.dat"
fname2 = "ref.dat"


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
fig.set_size_inches([15, 4])
plt.subplots_adjust(left=0.1, right=0.98, top=0.9, bottom=0.16)

# data read
data1=(np.loadtxt(fname1,dtype="float",skiprows=0)).T
print(np.shape(data1))

# data definition
jmax = 121
kmax = 201

x1=data1[0,:jmax]  # x
x2=data1[1,::jmax]  # y
q1=data1[2]  # r
q2=data1[3]  # u
q3=data1[4]  # v
q4=data1[5]  # w
q5=data1[6]  # p
q6=data1[7]  # ry1
# q6=data1[7]  # ry2
# q7=data1[8]  # ry3

q1 = np.array(q1).reshape(-1, jmax).tolist()
q2 = np.array(q2).reshape(-1, jmax).tolist()
q3 = np.array(q3).reshape(-1, jmax).tolist()
q4 = np.array(q4).reshape(-1, jmax).tolist()
q5 = np.array(q5).reshape(-1, jmax).tolist()
q6 = np.array(q6).reshape(-1, jmax).tolist()
# q7 = np.array(q7).reshape(-1, jmax).tolist()


#**  Blasius ********************
xat = 12
xpoint = int((xat/12)*(jmax-1))
uinf = 0.2

# eta = sym.symbols("eta")
# a = 0.332
# b = 1.73
# c = 0.231

# fdd = c*sym.exp(-0.25*(eta-b)**2)
# fd = [ -sym.integrate(fdd,(eta,e,INFINITY)) + 0.75 for e in x2/(xat-2)*np.sqrt(100)] 
# fd = -sym.integrate(fdd,(eta,eta,INFINITY))
# f = -sym.integrate(fd,(eta,eta,INFINITY))
# f = eta - b + f

fig.set_size_inches([10, 8])
data2=(np.loadtxt(fname2,dtype="float",skiprows=0)).T
print(np.shape(data2))

eta = data2[0]  # eta
y = data2[1]/15  # y
u = data2[2]  # u

# # *** 2D slice ****************************
plt.grid(False)
nmax = x2.shape[0]
plt.ylim(0,1.05)

plt.xlabel(r'$u/U$',fontsize=22,y=0)
plt.ylabel(r"$y$",fontsize=22)

u_uinf = np.array(q2).T[xpoint][:]/uinf
ref = np.ones(len(x2))

plt.plot(ref, x2, linewidth=1, color="gray", linestyle="dashed",zorder=1)
plt.plot(u_uinf, x2, linewidth=2, color="red", linestyle="solid",zorder=1)
plt.plot(u, y, linewidth=2, color="black", linestyle="solid",zorder=1)

ax = plt.gca()
aspect = 1/1* (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])                     
ax.set_aspect(aspect)

# *** 2D contour ***************************
# # plt.xlim(0,1)
# plt.ylim(0,0.4)
# plt.grid(False)
# nmax = x1.shape[0]

# # plt.title("Plot 2D array")
# plt.xlabel(r'$x$',fontsize=22)
# plt.ylabel(r"$y$",fontsize=22)
# # plt.xticks(np.arange(0,1.2,0.2))

# xx, yy = np.meshgrid(x1, x2)
# cf = plt.contourf(xx, yy, q1, cmap="jet")
# cf = plt.contourf(xx, yy, q2, cmap="jet")
# # cf = plt.contourf(xx, yy, q3, cmap="jet")
# # cf = plt.contourf(xx, yy, q4, cmap="jet")
# cf = plt.contourf(xx, yy, q5, cmap="jet",levels=36)
# # cf = plt.contourf(xx, yy, q6, cmap="jet")
# # plt.plot(x2,fd)

# cb = plt.colorbar(cf)
# cb.set_label(r'$u$',fontsize=22,y=1,x=1.2,rotation=0)
# cb.formatter.set_scientific(True)
# cb.ax.ticklabel_format(style='sci', scilimits=(-8,8)) 

# ax = plt.gca()
# aspect = 1/4* (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])                     
# ax.set_aspect(aspect)

# plt.show()
fig.savefig("2dtest.png")
