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
# fname1 = "output_200000.dat"
fname1 = "output_040000.dat"
# fname1 = "output_130000.dat"
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

q1 = np.array(q1).reshape(-1, jmax).tolist()
q2 = np.array(q2).reshape(-1, jmax).tolist()
q3 = np.array(q3).reshape(-1, jmax).tolist()
q4 = np.array(q4).reshape(-1, jmax).tolist()
q5 = np.array(q5).reshape(-1, jmax).tolist()
q6 = np.array(q6).reshape(-1, jmax).tolist()
# q7 = np.array(q7).reshape(-1, jmax).tolist()


#**  Blasius ********************
xat = 1  # observe point from wall start in x rigion
xpoint = int(((xat+2)/12)*(jmax-1))
uinf = 0.2
Rex  = 100 # @x=1

eta = x2 * np.sqrt(Rex/uinf)/np.sqrt(xat)  # eta = y/x * sq(Re*x)
u_one = np.ones(len(x2))
u_uinf = np.array(q2).T[xpoint][:]/uinf

d_thry = 5.3*np.sqrt(xat)/np.sqrt(Rex/uinf)
print("delta_thy : ",d_thry)

for i in range(len(x2)):
    if (u_uinf[i]>0.995) : 
        scale = i
        break

# scale = np.argmax(u_uinf)
d_calc = x2[scale+2]
print("delta_cal : ",d_calc)

#****** read reference solusion ************
data2=(np.loadtxt(fname2,dtype="float",skiprows=0)).T

d_ref = 5.32
eta2  = data2[0]   # eta
y     = data2[1]*(d_thry/d_ref) # fiten y with delta point
# y     = data2[1]/d_thry  # fiten y with delta point
ref_u = data2[2]  # u

x2 = x2/d_calc
# y  = y /d_calc

# x2 = x2/d_thry
y  = y /d_thry

# # *** 2D slice ****************************
fig.set_size_inches([10, 8])

plt.grid(False)
plt.xlabel(r'$u/U$',fontsize=22,y=0)
plt.ylabel(r"$y/\delta$",fontsize=22)
# plt.ylabel(r"$y$",fontsize=22)

plt.plot(u_one, x2, linewidth=1, color="gray", linestyle="dashed",zorder=1)

plt.ylim(0,2)
# plt.ylim(0,1)
plt.plot(u_uinf, x2, linewidth=2, color="red", linestyle="solid",zorder=1)
plt.plot(ref_u, y, linewidth=2, color="black", linestyle="solid",zorder=1)

# plt.ylim(0,8)
# plt.plot(u_uinf, eta , linewidth=2, color="red", linestyle="solid",zorder=1)
# plt.plot(ref_u , eta2, linewidth=2, color="black", linestyle="solid",zorder=1)

ax = plt.gca()
aspect = 1/1* (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])                     
ax.set_aspect(aspect)
fig.savefig("2dtest_udist.png")

# *** 2D contour ***************************
fig.set_size_inches([15, 4])
# plt.xlim(0,1)
# plt.ylim(0,10)
plt.grid(False)
nmax = x1.shape[0]

# plt.title("Plot 2D array")
plt.xlabel(r'$x$',fontsize=22)
# plt.ylabel(r"$y/\delta$",fontsize=22)
# plt.xticks(np.arange(0,1.2,0.2))

xx, yy = np.meshgrid(x1, x2)
cf = plt.contourf(xx, yy, q1, cmap="jet",levels=10)
cf = plt.contourf(xx, yy, q2, cmap="jet",levels=10)
# cf = plt.contourf(xx, yy, q3, cmap="jet")
# cf = plt.contourf(xx, yy, q4, cmap="jet")
# cf = plt.contourf(xx, yy, q5, cmap="jet",levels=36)
# cf = plt.contourf(xx, yy, q6, cmap="jet")
# plt.plot(x2,fd)
# plt.plot(x1, 5.3*np.sqrt(x1-2)/np.sqrt(Rex/uinf)/d_thry, linewidth=2, color="black", linestyle="solid",zorder=1)


cb = plt.colorbar(cf)
cb.set_label(r'$u$',fontsize=22,y=1,x=1.2,rotation=0)
cb.formatter.set_scientific(True)
cb.ax.ticklabel_format(style='sci', scilimits=(-8,8)) 

ax = plt.gca()
aspect = 1/4* (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])                     
ax.set_aspect(aspect)

# plt.show()
fig.savefig("2dtest.png")
