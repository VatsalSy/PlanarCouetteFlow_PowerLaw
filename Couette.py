# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

import numpy as np
import os
import subprocess as sp
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
matplotlib.rcParams['text.latex.unicode'] = True


def gettingfield(filename):
    print('Getting field values')
    exe = ["./getData", filename, str(xmin), str(xmax), str(ymin), str(ymax), str(nx), str(ny)]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    Xtemp = []
    Ytemp = []
    Utemp = []
    Vtemp = []
    D2temp = []
    for n1 in range(1,len(temp2)):
        temp3 = temp2[n1].split(" ")
        if temp3 == ['']:
            pass
        else:
            Xtemp.append(float(temp3[0]))
            Ytemp.append(float(temp3[1]))
            Utemp.append(float(temp3[2]))
            Vtemp.append(float(temp3[3]))
            D2temp.append(float(temp3[4]))
    X = np.asarray(Xtemp)
    Y = np.asarray(Ytemp)
    U = np.asarray(Utemp)
    V = np.asarray(Vtemp)
    D2 = np.asarray(D2temp)

    X.resize((nx+1, ny+1))
    Y.resize((nx+1, ny+1))
    U.resize((nx+1, ny+1))
    V.resize((nx+1, ny+1))
    D2.resize((nx+1, ny+1))

    print('Got field values')
    return X, Y, U, V, D2

# ----------------------------------------------------------------------------------------------------------------------


xmin = 0.0
xmax = 1.0
ymin = 0.0
ymax = 1.0
nx = 512
ny = 512

place = "lastBing"

name = "D2Bingham.png"
X, Y, U, V, D2 = gettingfield(place)
D2 /= np.sqrt(2.0)
## Part to plot
fig, ax = plt.subplots()
fig.set_size_inches(19.20, 10.80)
rc('axes', linewidth=2)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

# cntrl = ax.pcolormesh(X, Y, D2, cmap="jet",edgecolor='face')
cntrl = ax.pcolormesh(X, Y, D2, cmap="jet",edgecolor='face')
# cntrl = ax.pcolormesh(X, Y, V, cmap="jet",edgecolor='face')

cb1 = fig.add_axes([0.8, 0.1, 0.03, 0.8])
c1 = plt.colorbar(cntrl,cax=cb1)
# c1.set_label(r'$\|D_{ij}\|$', fontsize=30,labelpad=30)
c1.set_label(r'$U$', fontsize=30,labelpad=30)
# c1.set_label(r'$V$', fontsize=30,labelpad=30)
c1.ax.tick_params(labelsize=20)
c1.ax.tick_params(labelleft=True,labelright=False)

ax.set_xlabel(r'$X/D$', fontsize=30)
ax.set_ylabel(r'$Y/D$', fontsize=30)
ax.set_aspect('equal')
ax.set_ylim(xmin, xmax)
ax.set_xlim(ymin, ymax)

# plt.show()
plt.savefig(name,dpi=300)
plt.close()
