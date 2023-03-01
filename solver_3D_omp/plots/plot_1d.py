#!/bin/usr/python3
import numpy as np
import matplotlib.pyplot as plt

range_min=-1e-2
range_max= 1e-3
fname0 = "output_000000.dat"
fname1 = "output_001000.dat"

# global settings
plt.rcParams['xtick.direction'] = 'in'#x軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams['ytick.direction'] = 'in'#y軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams['xtick.major.width'] = 1.0 #x軸主目盛り線の線幅
plt.rcParams['ytick.major.width'] = 1.0 #y軸主目盛り線の線幅
plt.rcParams['font.size'] = 22 #フォントの大きさ
plt.rcParams['axes.linewidth'] = 1.0 # 軸の線幅edge linewidth。囲みの太さ
# plt.rcParams['text.usetex'] = True # Latex text
plt.rcParams["axes.formatter.use_mathtext"]=True

fig = plt.figure()
fig.set_size_inches([10, 8])
plt.subplots_adjust(left=0.12, right=0.95, top=0.9, bottom=0.13)

# data read
data0=(np.loadtxt(fname0,dtype="float",skiprows=0)).T
x0_1=data0[0]  # x
q0_1=data0[1]  # r
q0_2=data0[2]  # u
q0_3=data0[3]  # v
q0_4=data0[4]  # w
q0_5=data0[5]  # p
q0_6=data0[6]  # ry1
q0_7=data0[7]  # ry2
q0_8=data0[8]  # temp
q0_9=data0[9]*0.1  # gam


data1=(np.loadtxt(fname1,dtype="float",skiprows=0)).T
print(np.shape(data1))

# data definition
x1=data1[0]  # x
q1=data1[1]  # r
q2=data1[2]  # u
q3=data1[3]  # v
q4=data1[4]  # w
q5=data1[5]  # p
q6=data1[6]  # ry1
q7=data1[7]  # ry2
q8=data1[8]  # temp
q9=data1[9]*0.1  # gam

# format
plt.xlim(0,1)
plt.grid(False)
nmax = x1.shape[0]

plt.xlabel(r'$x$',fontsize=28,y=0)
plt.ylabel(r"$\rho Y_i,\: u, p, T$",fontsize=28)

plt.plot(x0_1[:nmax], q0_2[:nmax], "o", color="red"        , markevery=20, markersize="10", markerfacecolor="white", zorder=-2)
plt.plot(x0_1[:nmax], q0_5[:nmax], "o", color="green"      , markevery=20, markersize="10", markerfacecolor="white", zorder=-5)
plt.plot(x0_1[:nmax], q0_6[:nmax], "o", color="blue"       , markevery=20, markersize="10", markerfacecolor="white", zorder=-6)
plt.plot(x0_1[:nmax], q0_7[:nmax], "o", color="orange"     , markevery=20, markersize="10", markerfacecolor="white", zorder=-7)

plt.plot(x1[:nmax], q2[:nmax], linewidth=1.2, color="red"          , linestyle="solid",zorder=2)
plt.plot(x1[:nmax], q5[:nmax], linewidth=1.2, color="green"        , linestyle="solid",zorder=5)
plt.plot(x1[:nmax], q6[:nmax], linewidth=1.2, color="blue"         , linestyle="solid",zorder=6)
plt.plot(x1[:nmax], q7[:nmax], linewidth=1.2, color="orange"       , linestyle="solid",zorder=7)


# plt.show()
fig.savefig("1dtest.png")
# fig.savefig("1dtest.pdf")

