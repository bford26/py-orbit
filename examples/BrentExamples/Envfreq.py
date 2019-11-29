print "Start."
#-----------------------------------------------------------------------------
# IMPORTS
#-----------------------------------------------------------------------------

# Lattice Requirements
from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice, TEAPOT_MATRIX_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.matrix_lattice import MATRIX_Lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer

# Basic Utilities
from orbit.utils.consts import *
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit

# Distributions
from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import WaterBagDist2D, GaussDist3D, KVDist3D, SelfConDist2D, SelfConDist2D_2, SelfConDist3D

# TwissAnalysis
from orbit.diagnostics import diagnostics, diagnosticsLatticeModifications, TeapotDiagnosticsNode
from orbit.diagnostics.TeapotDiagnosticsNode import TeapotStatLatsNode

# SpaceCharge 2.5D
from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications
from spacecharge import SpaceChargeCalc2p5D, Boundary2D

# Beam Bunch
from bunch import Bunch


from Functions import BeamFunctions
from orbit.CustomNodes.CustomNodes import psDemonNode
from orbit.EnvelopeSolver.EnvelopeSolver import BasicSolverNode


import math
import random
import sys
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
#-----------------------------------------------------------------------------
# File Management
#-----------------------------------------------------------------------------
bf = BeamFunctions()
# Base folder where script runs
basePath = "/Users/6bf/py-orbit/examples/EnvelopeFrequency/"
# Here will be sent loop dat and the initial and final .dat files
datPath = basePath + "data/"
# This will move all plots to the folder
# plotPath = basePath + "plots/"
# This path will move the DemonNode .txt files here
txtPath = basePath + "dmNodes/"

bf.setDatPath(datPath)
# bf.setPlotPath(plotPath)
bf.setTextPath(txtPath)
bf.setBasePath(basePath)

#-----------------------------------------------------------------------------
# Test Options
#-----------------------------------------------------------------------------
EnvSolverNodes = True
SpaceCharge = True
DemonNodes = True

# This sets where Solver Child Nodes will be added
if EnvSolverNodes:
    Placement = AccNode.BODY

intensity = 1.e+14

INDEX = 1
mass = 0.93827231
macrosize = 1
KE = 1.0 # GeV kinetic energy

#How many lattice
# LatCount = 3

#Lengths of each lattice
Lattice_Length = 50.0
#How many Nodes in each lattice
NodeCount = 200

#Bunch count
Count = len(intensity)

#-----------------------------------------------------------------------------
# Twiss Parameters
#-----------------------------------------------------------------------------
alphaX = 0.0
alphaY = 0.0
betaX = Lattice_Length / ( 2*math.pi)
betaY = Lattice_Length / ( 2*math.pi)
emitX = 2.e-5
emitY = 2.e-5

#-----------------------------------------------------------------------------
# Envelope Bunch
#-----------------------------------------------------------------------------

def makeBunches(Count):
    Envs = []
    for i in range(Count):
        Env = Bunch()
        Env.mass(mass)
        Env.macroSize(macrosize)
        Env.getSyncParticle().kinEnergy(KE)

        a, ap = math.sqrt(betaX * emitX), 0.0
        b, bp = math.sqrt(betaY * emitY), 0.0
        SigE = 1.e-3

        Env.addParticle(a,ap,b,bp,0.,SigE)

        Envs.append(Env)
    return Envs

BunchEnvelopes = makeBunches(Count)

#-----------------------------------------------------------------------------
# Lattices
#-----------------------------------------------------------------------------
def AddUniformNode(lattice, elemNum, ElemLength, LenOfTune, x_tune, y_tune, LatType):
    """
    This Function will make a Uniform lattice node and
    add the created node to the lattice
    """
    if(elemNum < 10):
        index = "00" + str(elemNum)
    elif(elemNum < 100):
        index = "0" + str(elemNum)
    else:
        index = str(elemNum)

    elem = teapot.UniLatTEAPOT("UniTest" + "_Elem_" + index)
    lattice.addNode(elem)
    elem.setLength(ElemLength)
    elem.addParam("LenTunes", LenOfTune)
    elem.addParam("TuneX", x_tune)
    elem.addParam("TuneY", y_tune)
    elem.addParam("LatType", LatType)

LatticeType = "LINAC"
LenOfTune = Lattice_Length
xTune = 1.0
yTune = 1.0

if LatticeType == "LINAC":
    Dispersion = 0.0
    LType = 1
elif(LatticeType == "RING"):
    Dispersion = betaX**2 * 2*pi / Lattice_Length
    LType = 0

# LatList = []

# for i in range(Count):
    # if(i < 100 and i < 10):
    #     tag = "00"+str(i)
    # elif(i < 100):
    #     tag = "0"+str(i)
    # else:
    #     tag = str(i)

lat = teapot.TEAPOT_Lattice("Lattice_")# + tag)
    # LatList.append(lat)

i=1
# for i,lat in enumerate(LatList):
for j in range(NodeCount):
    AddUniformNode(lat, j, Lattice_Length/NodeCount, LenOfTune, xTune, yTune, LType)


# for i in range(LatCount):

if EnvSolverNodes:
    ChildNode = BasicSolverNode()
    ChildNode.setDispersion(Dispersion)
    ChildNode.setDLength(0.0001)
    for node in lat.getNodes():
        node.addChildNode(ChildNode, Placement)

if DemonNodes:

    dmlist = []

    for i,node in enumerate(lat.getNodes()):
        if(i < 10 and i < 100):
            index = "00" + str(i)
        elif(i < 100):
            index = "0" + str(i)
        else:
            index = str(i)

        DemonNodePlacement = AccNode.EXIT
        FileName = "dm" + index
        dmlist.append(FileName)
        node.addChildNode(psDemonNode(FileName,fdump=True, lloop=True), DemonNodePlacement)

# File Management
for i,file in enumerate(dmlist):
    dmlist[i] = file + ".txt"





#Frequency Function
E = 4 * emitX
#Calculate spacecharge

OMEGA = []
SPACECHARGE = []

bunch = BunchEnvelopes[INDEX]

beta = bunch.getSyncParticle().beta()
gamma = bunch.getSyncParticle().gamma()
r_p = proton_classical_radius
# beam partciles/ length
lam = intensity[INDEX] / Lattice_Length
SC = 2 * r_p * lam / ( beta**2 * gamma**3)
# SPACECHARGE.append(SC)

paramsDict = {}
paramsDict["bunch"] = bunch
paramsDict["EmitX"] = emitX
paramsDict["EmitY"] = emitY
paramsDict["spacecharge"] = SC


lat.trackBunch(bunch, paramsDict)

bf.MoveFile(dmlist, bf.TXTpath)


#
# for i,file in enumerate(dmlist):
#     data = bf.ReadData(file)
#     x[i] = data[0]
#     y[i] = data[2]


# INDEX = 138
# OMEGA = 2 * math.pi / ( INDEX / 200.0 * Lattice_Length )
#
# E = 4 * emitX
# TuneTerm = (2 * pi * xTune / Lattice_Length)**2
# scline = np.linspace(0, SC*1.05, 2000)
# omega = 2 * np.sqrt( TuneTerm - 0.5*scline/(2*E)*np.sqrt( TuneTerm + (scline*scline)/(4*E*E)) + (scline*scline)/(4*E*E) )

# plt.figure()
# plt.plot(scline,omega)
# plt.scatter([SC], [OMEGA], c='r')
#
#
# plt.savefig("xplots.png")



#
# for j,bunch in enumerate(BunchEnvelopes):
#
#     beta = bunch.getSyncParticle().beta()
#     gamma = bunch.getSyncParticle().gamma()
#     r_p = proton_classical_radius
#     # beam partciles/ length
#     lam = intensity[j] / Lattice_Length
#     SC = 2 * r_p * lam / ( beta**2 * gamma**3)
#
#     SPACECHARGE.append(SC)
#
#
#     paramsDict = {}
#     paramsDict["bunch"] = bunch
#     paramsDict["EmitX"] = emitX
#     paramsDict["EmitY"] = emitY
#     paramsDict["spacecharge"] = SC
#
#     lat.trackBunch(bunch, paramsDict)
#
#     bf.MoveFile(dmlist, bf.TXTpath)
#
#     x = np.zeros((NodeCount,1))
#     y = np.zeros((NodeCount,1))
#
    # for i,file in enumerate(dmlist):
    #     data = bf.ReadData(file)
    #     x[i] = data[0]
    #     y[i] = data[2]
#
#     print(x)
#     plt.subplot(len(BunchEnvelopes), 1, j+1)
#     plt.plot(x)
#
#     plt.savefig("xplots.png")
#     plt.close()




    # xmax = 0
    # for i in range(len(x)):
    #     if xmax < x[i]:
    #         xmax = x[i]
    #         index = i
    #
    #
    # for i in range(len(x)):
    #     if(x[i] >= xmax*(0.999)  and i != index):
    #         matchIndex = i
    #         # print("====="+str(j)+"=====")
    #         # print(matchIndex)
    #         # print(x[matchIndex])
    #
    #         if j == 3:
    #             S = 137 * Lattice_Length / NodeCount
    #
    #             PERIOD = S
    #             W_env = 2 * pi / PERIOD
    #             OMEGA.append(W_env)





# plt.title("Frequencies")
#
# #THIS will plot the prediction
#
#
#
#
#

#
#
# #plot the data from simulation
# # for i in range(Count):
#
#
# plt.xlabel("SPACECHARGE")
# plt.ylabel("OMEGA")
# #
# plt.savefig("FrequencyData.png")
# plt.close()
#

print "Stop."

#EOF
