#!/bin/usr/python3
import numpy as np
import matplotlib.pyplot as plt

range_min=-1e-2
range_max= 1e-3
fname0 = "norm_div.dat"
fname1 = "norm_pro.dat"
fname2 = "norm_pro2.dat"
fname3 = "norm_keeppe.dat"

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
plt.subplots_adjust(left=0.12, right=0.95, top=0.9, bottom=0.13)

# data read
data0=(np.loadtxt(fname0,dtype="float",skiprows=0)).T
data1=(np.loadtxt(fname1,dtype="float",skiprows=0)).T
data2=(np.loadtxt(fname2,dtype="float",skiprows=0)).T
data3=(np.loadtxt(fname3,dtype="float",skiprows=0)).T

# data definition
x0 = data0[0]/10
x1 = data1[0]/10
x2 = data2[0]/10
x3 = data3[0]/10
i = 7
q0 = data0[i]
q1 = data1[i]
q2 = data2[i]
q3 = data3[i]
q5=np.zeros(len(x1))  # exact

# format
# plt.xlim(0,1)
plt.grid(False)
plt.yscale("log")

plt.xlabel(r'$t$',fontsize=28,y=0)
plt.ylabel(r"$L2(p)$",fontsize=28)

plt.plot(x0, q0, linewidth=1.2, color="blue"  , linestyle="solid",zorder=1)
# plt.plot(x1, q1, linewidth=1.2, color="green"  , linestyle="solid",zorder=1)
plt.plot(x2, q2, linewidth=1.2, color="green"   , linestyle="solid",zorder=2)
plt.plot(x3, q3, linewidth=1.2, color="red"  , linestyle="solid",zorder=3)


# plt.show()
fig.savefig("step.png")
# fig.savefig("1dtest.pdf")

