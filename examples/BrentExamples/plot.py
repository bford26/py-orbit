print "Start."
#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from mpl_toolkits import mplot3d
import matplotlib.figure

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

from BFunctions import BeamFunctions

from bunch import Bunch

#-----------------------------------------------------------------------------
# File Management
#-----------------------------------------------------------------------------
bf = BeamFunctions()

bf.setTextPath("/Users/6bf/py-orbit/examples/BrentExample/dmNodes/")
bf.setBasePath("/Users/6bf/py-orbit/examples/BrentExample/")
bf.setPlotPath("/Users/6bf/py-orbit/examples/BrentExample/plots/")
bf.setDatPath("/Users/6bf/py-orbit/examples/BrentExample/data/")

bf.setSpaceCharge(True) # No spacecharge
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
# filelist = bf.FindFiles(".txt")
#
# for file in filelist:
#     bf.plotANY(file)


data = bf.ReadData("initial.dat")

f = plt.figure(dpi=600)
ax1 = f.add_subplot(111)

ax1.set_title("Initial Distribution")
ax1.set_xlabel("x [mm]", fontsize=14)
ax1.set_ylabel("y [mm]", fontsize=14)


x = 1000*data[0]
y = 1000*data[2]

ax1.scatter(x,y, s=1)

xline = [ np.linspace(0.0,max(x),200), np.zeros((200,)) ]
yline = [ np.zeros((200,)), np.linspace(0.0,max(y),200) ]

ax1.plot(xline[0],xline[1], 'k', linewidth=2.0)
ax1.plot(yline[0],yline[1], 'k', linewidth=2.0)

ax1.set_xlim([-1.25*max(x),1.25*max(x)])
ax1.set_ylim([-1.3*max(y), 1.3*max(y)])

txta = "a"
txtb = "b"
ax1.annotate(txta, xy=(max(x), 0.0), xytext=(max(x)*1.11,0.0), horizontalalignment='right', verticalalignment='bottom', fontsize=14)
ax1.annotate(txtb, (0.0, max(y)), xytext=(0.0,max(y)*1.2), horizontalalignment='right', verticalalignment='top', fontsize=14)

f.savefig("envelope_plot.png", dpi=600)
plt.close()


#
#
# f = open("BetaData_nosc.txt")
#
# x = []
# y = []
#
# for line in f.readlines():
#     tokens = line.split()
#     x.append(float(tokens[0]))
#     y.append(float(tokens[1]))
#
# f.close()
#
# nosc = [x, y]
#
# f = open("BetaData_sc.txt")
#
# x = []
# y = []
#
# for line in f.readlines():
#     tokens = line.split()
#     x.append(float(tokens[0]))
#     y.append(float(tokens[1]))
#
# f.close()
#
# sc = [x, y]
#
# f = plt.figure(dpi=600)
#
#
# ax1 = f.add_subplot(211, ylabel="[m/rad]")#, ylim=(0.99*min(nosc[1]),1.01*max(nosc[1])))
# ax2 = f.add_subplot(212, ylabel="[m/rad]", xlabel="Node")#, ylim=(0.99*min(sc[1]),1.01*max(sc[1])))
#
# ax1.set_title("Statistical Beta Functions")
#
# ax1.plot(nosc[0], label="Beta X")
# ax1.plot(nosc[1], label="Beta Y")
#
# txt1 = "No SpaceCharge"
# ax1.text(0.05,0.95, txt1, transform=ax1.transAxes, verticalalignment='top', bbox=dict(boxstyle="square",facecolor="white"))
# ax1.legend(loc="upper right")
#
# ax2.plot(sc[0],label="Beta X")
# ax2.plot(sc[1],label="Beta Y")
# ax2.legend(loc="upper right")
#
#
# txt2 = "SpaceCharge"
# ax2.text(0.05,0.95, txt2, transform=ax2.transAxes, verticalalignment='top', bbox=dict(boxstyle="square",facecolor="white"))
#
# f.savefig("BetaPlots.png", dpi=600)
# plt.close()






# bf.betaplots()




# bf.subplots("initial.dat","final.dat", "initvsfin.png", "Distributions",scatter=True)

# bf.subplots("betaData_nosc.txt", "betaData_sc.txt", "betaplots.png", "Beta Function" )

#
# ftop = "betaData_nosc.txt"
# fbot = "betaData_sc.txt"
#
# fpath = self.basepath
# f = open(fpath+ftop)
#
# x = []
# x = []
#
# for line in f.readlines():
#     token = line.split()
#     x1 = float(token[0])
#     y1 = float(token[1])
#     x.append(x1)
#     y.append(y1)
#
#
#































#-----------------------------------------------------------------------------
# Basic Phase Plot
#-----------------------------------------------------------------------------

# bf.plotANY("")



# bf.clearPlots()
#-----------------------------------------------------------------------------
# X-Y phase plane for all psDemonNode Files
#-----------------------------------------------------------------------------

# filelist = bf.FindFiles(".txt")
# for file in filelist:
#     bf.plotANY(file,"xy")

# bf.plotANY("initial.dat","xy")


#-----------------------------------------------------------------------------
# Inital and Final All Phase Plots
#-----------------------------------------------------------------------------

# phases = ["xy"]#, "xxp", "yyp", "xyp", "yxp"]
# for phase in phases:
#     bf.plotANY("initial.dat", phase)
#     bf.plotANY("final.dat", phase)

#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Single Particle plot through lattice
#-----------------------------------------------------------------------------
#
# filelist = bf.FindFiles(".txt")
#
# nParticles = 1
# nElements = len(filelist)
#
# x = np.zeros((nElements,nParticles))
# xp = np.zeros((nElements,nParticles))
# y = np.zeros((nElements,nParticles))
# yp = np.zeros((nElements,nParticles))
#
# i=0
# for file in filelist:
#
#     data = bf.ReadData(bf.TXTpath + file)
#
#     x[i] = float(data[0])
#     xp[i] = float(data[1])
#     y[i] = float(data[2])
#     yp[i] = float(data[3])
#     i += 1
#
# plt.figure()
# plt.title("x, y - plot")
# plt.subplot(2,1,1)
# plt.plot(x)
#
# plt.subplot(2,1,2)
# plt.plot(y)
#
# plt.savefig("xandyplot.png")
#




print "Stop."


















#eof
