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

from orbit.EnvelopeSolver.EnvelopeSolver import BasicSolverNode, TiltSolverNode


#-----------------------------------------------------------------------------
# File Pathing
#-----------------------------------------------------------------------------
bf = BeamFunctions()

# Base folder where script runs
basePath = "/Users/6bf/py-orbit/examples/TiltingEnvelope/"
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

RunParticleBunch = False


# This sets where Solver Child Nodes will be added
Placement = AccNode.BODY
DemonNodePlacement = AccNode.EXIT



Lattice_Length = 50.0
NodeCount = 200
intensity = 1.0e+14

mass = 0.93827231
macrosize =  intensity
KE = 1.0 # GeV kinetic energ

nParticles = 10000


def printLattice(lat):
    i=0
    for node in lat.getNodes():
        print("Parent Node: "+ node.getType())

        for cnode in node.getChildNodes(AccNode.ENTRANCE):
            print("\tEntrance: " + cnode.getType())

        for cnode in node.getBodyChildren():
            print("\tBody: " + cnode.getType())

        for cnode in node.getChildNodes(AccNode.EXIT):
            print("\tExit: " + cnode.getType())

        if(i>10):
            break
        i += 1


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
Env.macroSize(intensity/2)
Env.getSyncParticle().kinEnergy(KE)


base_radius = 2 * math.sqrt(betaX * emitX)

a = base_radius*1
e = base_radius*0.25
ap, ep = e / math.sqrt(betaX * betaY)  , -1.0 * a / math.sqrt(betaX * betaY)

b = base_radius*0.5
f = base_radius*1
bp, fp = f / math.sqrt(betaX * betaY) , -1.0 * b / math.sqrt(betaX * betaY)

Env.addParticle(a,ap,e,ep,0.0,0.0)
Env.addParticle(b,bp,f,fp,0.0,0.0)

#-----------------------------------------------------------------------------
# Self-con Bunch
#-----------------------------------------------------------------------------
TwissX = TwissContainer(alphaX, betaX, emitX)
TwissY = TwissContainer(alphaY, betaY, emitY)


if RunParticleBunch:
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

if RunParticleBunch:
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
Xtune = 1.0
Ytune = 1.0
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
    if RunParticleBunch:
        MakeUniformNodes(BLat, i, ElementLength, LengthOfTune, Xtune, Ytune, LType)

#-----------------------------------------------------------------------------
# Child Nodes
#-----------------------------------------------------------------------------
if EnvSolverNodes:
    dL = Lattice_Length / NodeCount
    for node in EnvLat.getNodes():
        node.addChildNode(TiltSolverNode(length=dL), Placement)

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

    if RunParticleBunch:
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
if SpaceCharge == True and RunParticleBunch == True:
    sizeX = 64  #number of grid points in horizontal direction
    sizeY = 64  #number of grid points in vertical direction
    sizeZ = 32  #number of longitudinal slices in the 3D space charge solver

    calc2p5d = SpaceChargeCalc2p5D(sizeX, sizeY, sizeZ)
    sc_path_length_min = 0.00001 # Not 100% on this particular value
    scLatticeModifications.setSC2p5DAccNodes(BLat, sc_path_length_min, calc2p5d)

# if RunParticleBunch == True:

    # doesnt do the tilting thing i want
    # tilting_angle = np.linspace(0,2*pi,NodeCount*2)
    # i=0
    # for node in BLat.getNodes():
    #
    #     for cnode in node.getChildNodes(AccNode.ENTRANCE):
    #         if cnode.getType() == "tilt teapot":
    #             cnode.setTiltAngle(tilting_angle[i])
    #             i+=1
    #     for cnode in node.getChildNodes(AccNode.EXIT):
    #         if cnode.getType() == "tilt teapot":
    #             cnode.setTiltAngle(tilting_angle[i])
    #             i+=1


#-----------------------------------------------------------------------------
# Track
#-----------------------------------------------------------------------------
Env.dumpBunch("initial_Env.dat")

if RunParticleBunch:
    b.dumpBunch("initial_B.dat")

    # beta = b.getSyncParticle().beta()
    # gamma = b.getSyncParticle().gamma()
    # r_p = proton_classical_radius
    # b_SC = 2 * r_p * lam / (beta**2 * gamma**3)

#Calculate spacecharge
beta = Env.getSyncParticle().beta()
gamma = Env.getSyncParticle().gamma()
r_p = proton_classical_radius
# beam partciles / length = lineary charge density
lam = intensity / Lattice_Length
Env_SC = 2 * r_p * lam / (beta**2 * gamma**3)



paramsDict = {}
paramsDict["bunch"] = Env
paramsDict["spacecharge"] = Env_SC
paramsDict["dlength"] = dL
# params

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



if RunParticleBunch:
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

if RunParticleBunch:
    b.dumpBunch("final_B.dat")



#-----------------------------------------------------------------------------
# File Management
#-----------------------------------------------------------------------------
if RunParticleBunch:
    files = ["initial_B.dat","final_B.dat"]
    bf.MoveFile(files, bf.DATpath)

files = ["initial_Env.dat","final_Env.dat"]
bf.MoveFile(files, bf.DATpath)


# for node in BLat.getNodes():
#     i=0
#     for cnode in node.getChildNodes(AccNode.ENTRANCE):
#         if cnode.getType() == "tilt teapot":
#             print("In: "+str(cnode.getTiltAngle()))
#
#     for cnode in node.getChildNodes(AccNode.EXIT):
#         if cnode.getType() == "tilt teapot":
#             print("Out: "+str(cnode.getTiltAngle()))
#
#     if i>10:
#         break
#     i+=1
#


# printLattice(EnvLat)

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
