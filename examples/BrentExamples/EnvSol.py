print "Start."
#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import math
import random
import sys
import os

import numpy as np
import matplotlib.pyplot as plt

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
#-----------------------------------------------------------------------------
from Functions import BeamFunctions, EnvelopeFunctions
from orbit.CustomNodes.CustomNodes import psDemonNode

from orbit.EnvelopeSolver.EnvelopeSolver import BasicSolverNode


#-----------------------------------------------------------------------------
# File Pathing
#-----------------------------------------------------------------------------
bf = BeamFunctions()

# Base folder where script runs
basePath = "/Users/6bf/py-orbit/examples/EnvelopeSolver/"
# Here will be sent loop dat and the initial and final .dat files
datPath = basePath + "data/"
# This will move all plots to the folder
plotPath = basePath + "plots/"
# This path will move the DemonNode .txt files here
txtPath = basePath + "dmNodes/"

bf.setDatPath(datPath)
bf.setPlotPath(plotPath)
bf.setTextPath(txtPath)
bf.setBasePath(basePath)

# bf.RemoveAllData()

#-----------------------------------------------------------------------------
# Test Options
#-----------------------------------------------------------------------------
EnvSolverNodes = True
DemonNodes = True
SpaceCharge = True

# This sets where Solver Child Nodes will be added
Placement = AccNode.BODY
DemonNodePlacement = AccNode.EXIT

Lattice_Length = 50.0
NodeCount = 200
intensity = 1.0e+14

mass = 0.93827231
macrosize =  intensity
KE = 1.0 # GeV kinetic energ

nParticles = 100000

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
Env = Bunch()
Env.mass(mass)
Env.macroSize(macrosize)
Env.getSyncParticle().kinEnergy(KE)


base_radius = 2 * math.sqrt(betaX * emitX)

a = base_radius
b = a
ap, bp = 0.0 , 0.0
SigE = 1.e-3

Env.addParticle(a,ap,b,bp,0.,SigE)


#-----------------------------------------------------------------------------
# Self-con Bunch
#-----------------------------------------------------------------------------
TwissX = TwissContainer(alphaX, betaX, emitX)
TwissY = TwissContainer(alphaY, betaY, emitY)

distro = SelfConDist2D_2(TwissX, TwissY)

b = Bunch()
b.mass(mass)
b.macroSize(intensity/nParticles)
b.getSyncParticle().kinEnergy(KE)

for i in range(nParticles):
    (x,xp,y,yp) = distro.getCoordinates()
    z = (2*random.random() - 1) * Lattice_Length
    zp = 0.0
    b.addParticle(x,xp,y,yp,z,zp)



#-----------------------------------------------------------------------------
# Create Lattice
#-----------------------------------------------------------------------------
EnvLat = teapot.TEAPOT_Lattice("Envelope Solition Lattice")
BLat = teapot.TEAPOT_Lattice("Bunch Test Lattice")

def MakeUniformNodes(TLattice, elemNum, ElemLength, LenOfTune, x_tune, y_tune, LatType):
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
    TLattice.addNode(elem)
    elem.setLength(ElemLength)
    elem.addParam("LenTunes", LenOfTune)
    elem.addParam("TuneX", x_tune)
    elem.addParam("TuneY", y_tune)
    elem.addParam("LatType", LatType)

#-----------------------------------------------------------------------------
# Uniform Lattice Node Parameters
#-----------------------------------------------------------------------------
ElementLength = Lattice_Length / NodeCount
LengthOfTune = Lattice_Length
Xtune = 1.
Ytune = 1.
LatticeType = "LINAC"
# LatticeType = "RING"

if LatticeType == "LINAC":
    Dispersion = 0.0
    LType = 1
elif(LatticeType == "RING"):
    Dispersion = betaX**2 * 2*pi / Lattice_Length
    LType = 0

for i in range(NodeCount):
    MakeUniformNodes(EnvLat, i, ElementLength, LengthOfTune, Xtune, Ytune, LType)
    MakeUniformNodes(BLat, i, ElementLength, LengthOfTune, Xtune, Ytune, LType)

#-----------------------------------------------------------------------------
# Child Nodes
#-----------------------------------------------------------------------------
if EnvSolverNodes:
    dL = Lattice_Length / NodeCount
    Disp = Dispersion
    for node in EnvLat.getNodes():
        node.addChildNode(BasicSolverNode(dL, Disp), Placement)

if DemonNodes:
    i = 0
    for node in EnvLat.getNodes():
        if(i < 10 and i < 100):
            index = "00" + str(i)
        elif(i < 100):
            index = "0" + str(i)
        else:
            index = str(i)

        FileName = "dm" + index
        node.addChildNode(psDemonNode(FileName,fdump=True), DemonNodePlacement)
        i += 1

    i = 0
    for node in BLat.getNodes():
        if(i < 10 and i < 100):
            index = "00" + str(i)
        elif(i < 100):
            index = "0" + str(i)
        else:
            index = str(i)

        FileName = "dm" + index
        node.addChildNode(psDemonNode(FileName,fdump=True), DemonNodePlacement)
        i += 1

#-----------------------------------------------------------------------------
# SpaceCharge for Bunch Lattice
#-----------------------------------------------------------------------------
if SpaceCharge == True:
    sizeX = 64  #number of grid points in horizontal direction
    sizeY = 64  #number of grid points in vertical direction
    sizeZ = 32  #number of longitudinal slices in the 3D space charge solver

    calc2p5d = SpaceChargeCalc2p5D(sizeX, sizeY, sizeZ)
    sc_path_length_min = 0.00001 # Not 100% on this particular value
    scLatticeModifications.setSC2p5DAccNodes(BLat, sc_path_length_min, calc2p5d)

#-----------------------------------------------------------------------------
# Track
#-----------------------------------------------------------------------------
Env.dumpBunch("initial_Env.dat")
b.dumpBunch("initial_B.dat")

#Calculate spacecharge
beta = Env.getSyncParticle().beta()
gamma = Env.getSyncParticle().gamma()
r_p = proton_classical_radius
# beam partciles / length = lineary charge density
lam = intensity / Lattice_Length
Env_SC = 2 * r_p * lam / (beta**2 * gamma**3)

beta = b.getSyncParticle().beta()
gamma = b.getSyncParticle().gamma()
b_SC = 2 * r_p * lam / (beta**2 * gamma**3)

paramsDict = {}
paramsDict["bunch"] = Env
paramsDict["EmitX"] = emitX
paramsDict["EmitY"] = emitY
paramsDict["spacecharge"] = Env_SC
paramsDict["dlength"] = dL
paramsDict["dispersion"] = Dispersion

EnvLat.trackBunch(Env, paramsDict)

# Moving Files
if DemonNodes == True:
    for i in range(NodeCount):
        if(i < 10 and i < 100):
            index = "00" + str(i)
        elif(i < 100):
            index = "0" + str(i)
        else:
            index = str(i)
        name = "dm" + index + ".txt"
        bf.MoveFile(name, bf.TXTpath + "EnvLat/")



paramsDict["bunch"] = b
BLat.trackBunch(b, paramsDict)

# Moving Files
if DemonNodes == True:
    for i in range(NodeCount):
        if(i < 10 and i < 100):
            index = "00" + str(i)
        elif(i < 100):
            index = "0" + str(i)
        else:
            index = str(i)
        name = "dm" + index + ".txt"
        bf.MoveFile(name, bf.TXTpath + "BLat/")


Env.dumpBunch("final_Env.dat")
b.dumpBunch("final_B.dat")



#-----------------------------------------------------------------------------
# File Management
#-----------------------------------------------------------------------------
files = ["initial_B.dat","final_B.dat"]
bf.MoveFile(files, bf.DATpath)

files = ["initial_Env.dat","final_Env.dat"]
bf.MoveFile(files, bf.DATpath)



# ef = EnvelopeFunctions()
# basePath = "/Users/6bf/py-orbit/examples/EnvelopeSolver/"
# # Here will be sent loop dat and the initial and final .dat files
# datPath = basePath + "data/"
# # This will move all plots to the folder
# plotPath = basePath + "plots/"
# # This path will move the DemonNode .txt files here
# txtPath = basePath + "dmNodes/"
# ef.setDatPath(datPath)
# ef.setPlotPath(plotPath)
# ef.setTextPath(txtPath)
# ef.setBasePath(basePath)


# ef.PlotAmplitudes()


print "Stop."
#eof
