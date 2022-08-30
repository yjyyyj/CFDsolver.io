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

fname1 = "output_000000.dat"
files = np.sort(list(pathlib.Path('./').glob('output_*.dat')))
fnames = [ i.name for i in files]
print(fnames)

# global settings
plt.rcParams['xtick.direction'] = 'in'#x軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams['ytick.direction'] = 'in'#y軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams['xtick.major.width'] = 1.0 #x軸主目盛り線の線幅
plt.rcParams['ytick.major.width'] = 1.0 #y軸主目盛り線の線幅
plt.rcParams['font.size'] = 12 #フォントの大きさ
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

    q1 = np.array(q1).reshape(-1, jmax)
    q2 = np.array(q2).reshape(-1, jmax)
    q3 = np.array(q3).reshape(-1, jmax)
    q4 = np.array(q4).reshape(-1, jmax)
    q5 = np.array(q5).reshape(-1, jmax)
    q6 = np.array(q6).reshape(-1, jmax)
    q7 = np.array(q7).reshape(-1, jmax)
    
    # format
    nmax = x1.shape[0]

    xx, yy = np.meshgrid(x1, x2)
    ax = fig.add_subplot(111, projection='3d')

    # ax.set_xlim(0,1)
    # ax.set_ylim(range_min,range_max)
    # ax.set_zlim(range_min,range_max)

    # plt.title("Plot 2D array")
    ax.set_xlabel(r'$x$',fontsize=28,y=0)
    ax.set_ylabel(r"$y$",fontsize=28)

    # ax.set_xticks([])
    # ax.set_yticks([])
    ax.set_zticks([])
    ax.zaxis.line.set_color((0.,0.,0.,0.))
    ax.w_xaxis.set_pane_color((0.,0.,0.,0.))
    ax.w_yaxis.set_pane_color((0.,0.,0.,0.))
    ax.w_zaxis.set_pane_color((0.,0.,0.,0.))
    ax.grid(False)

    surf = ax.plot_surface(xx, yy, q1, cmap='jet',rstride=1,cstride=1,linewidth=0,antialiased=False, shade=True)
    # surf = ax.plot_surface(xx, yy, q2, cmap='jet',rstride=1,cstride=1,linewidth=0,antialiased=False, shade=True)
    # surf = ax.plot_surface(xx, yy, q3, cmap='jet',rstride=1,cstride=1,linewidth=0,antialiased=False, shade=True)
    # surf = ax.plot_surface(xx, yy, q4, cmap='jet',norm=Normalize(vmin=range_min, vmax=range_max),rstride=1,cstride=1,linewidth=0,antialiased=False, shade=True)
    # surf = ax.plot_surface(xx, yy, q5, cmap='jet',rstride=1,cstride=1,linewidth=0,antialiased=False, shade=True)
    # surf = ax.plot_surface(xx, yy, q6, cmap='jet',rstride=1,cstride=1,linewidth=0,antialiased=False, shade=True)

    cb = fig.colorbar(surf)
    cb.set_label(r'$\rho$',fontsize=28,y=1,x=1.2,rotation=0)
    cb.formatter.set_scientific(True)
    # ax.view_init(elev=60, azim=45)
    cb.ax.ticklabel_format(style='sci', scilimits=(-3,3)) 

    # plt.show()
    fig.savefig("img{:03d}.png".format(i))

fp.run(stream)