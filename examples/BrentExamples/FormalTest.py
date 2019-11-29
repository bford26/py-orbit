print "Start."
#-----------------------------------------------------------------------------
# imports
#-----------------------------------------------------------------------------
import math
import random
import sys
import os
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
## Use for tracking thorughout the lattice
from orbit.CustomNodes.CustomNodes import psDemonNode
## Extra functions to handle reading data and exporting for plots
from BFunctions import BeamFunctions

#-----------------------------------------------------------------------------
# File Management
#-----------------------------------------------------------------------------
# Here needs to be basic file path the folder is located
bf = BeamFunctions()

# Base folder where script runs
basePath = "/Users/6bf/py-orbit/examples/BrentExample/"
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

bf.RemoveAllData()

#-----------------------------------------------------------------------------
# Data Printing Options
#-----------------------------------------------------------------------------
PrintLoop = False #Makes loop.dat files
PrintLatticeNodes = False #Prints all Nodes and ChildNodes for Lattice
NodeTrack = True #Will dump the bunch attributes at each node in file
PlotBeta = False


#-----------------------------------------------------------------------------
# Simulation Options
#-----------------------------------------------------------------------------
SpaceCharge = True #Include SpaceCharge
DemonNodes = True #Used for tracking nodal information
RunLoops = False #Actually loop for specified amount otherwise only 2 loops

SingleParticle = False #When active it puts only one particle in the lattice

#-----------------------------------------------------------------------------

nloops = 5
Lattice_Length = 50.0 #Meters
NodeCount = 200 #Lattice Elements/Nodes

nParticles = 10000 #Number of particles in the bunch
total_macroSize = 1.e+14 #Actual Number of particles being simulated


if SingleParticle ==True:
    nParticles = 1
    macrosize = nParticles
    PlotBeta = False
else:
    macrosize = total_macroSize / nParticles

#-----------------------------------------------------------------------------
# Defining Bunch Attributes
#-----------------------------------------------------------------------------
energy = 1.0 #GeV
b = Bunch()
b.mass(0.93827231)
b.macroSize(macrosize)
b.getSyncParticle().kinEnergy(energy)

#-----------------------------------------------------------------------------
# Twiss Parameters for Distribution
#-----------------------------------------------------------------------------
alphaX = 0.0
alphaY = 0.0
betaX = Lattice_Length / (2 * pi)
betaY = Lattice_Length / (2 * pi)
emittanceX = 2.e-5
emittanceY = 2.e-5

twissX = TwissContainer(alphaX, betaX, emittanceX)
twissY = TwissContainer(alphaY, betaY, emittanceY)

#-----------------------------------------------------------------------------
# Making Distribution
#-----------------------------------------------------------------------------
# dist = SelfConDist2D(twissX, twissY) #Code wrote by Austin H.
dist = SelfConDist2D_2(twissX, twissY) #Code wrote by J. Holmes adapted to python

if SingleParticle==True:
    PlotBeta = False
    x, xp = 0.0795, 0.
    y, yp = 0., -0.01
    z, zp = 0., 0.
    b.addParticle(x,xp,y,yp,z,zp)

else:
    for i in range(nParticles):
        (x,xp,y,yp) = dist.getCoordinates()
        z = ( 2 * random.random() -  1 ) * ( 0.5 * Lattice_Length )
        zp = 0.
        b.addParticle(x,xp,y,yp,z,zp)

#Print bunch data to file
b.dumpBunch("initial.dat")

paramsDict = {}
paramsDict["bunch"] = b
# lostbunch = Bunch()
# lostbunch.addPartAttr("LostParticleAttributes")
# paramsDict["lostbunch"] = lostbunch

#-----------------------------------------------------------------------------
# Create Lattice and Nodes
#-----------------------------------------------------------------------------
Lattice = teapot.TEAPOT_Lattice("Uniform Lattice Test")

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

#-----------------------------------------------------------------------------

if LatticeType == "RING":
    LType = 0
elif LatticeType == "LINAC":
    LType = 1

for i in range(NodeCount):
    MakeUniformNodes(Lattice, i, ElementLength, LengthOfTune, Xtune, Ytune, LType)


#-----------------------------------------------------------------------------
# Add SpaceCharge
#-----------------------------------------------------------------------------
bf.setSpaceCharge(SpaceCharge)
if SpaceCharge == True:
    sizeX = 64  #number of grid points in horizontal direction
    sizeY = 64  #number of grid points in vertical direction
    sizeZ = 32  #number of longitudinal slices in the 3D space charge solver

    calc2p5d = SpaceChargeCalc2p5D(sizeX, sizeY, sizeZ)
    sc_path_length_min = 0.000001 # Not 100% on this particular value
    scLatticeModifications.setSC2p5DAccNodes(Lattice, sc_path_length_min, calc2p5d)

#-----------------------------------------------------------------------------
# Add psDemonNodes for tracking bunch at the EXIT of each node
#-----------------------------------------------------------------------------
if DemonNodes == True:
    i = 0
    for node in Lattice.getNodes():
        if(i < 10 and i < 100):
            index = "00" + str(i)
        elif(i < 100):
            index = "0" + str(i)
        else:
            index = str(i)

        DemonNodePlacement = AccNode.BODY
        FileName = "dm" + index
        node.addChildNode(psDemonNode(FileName, fdump=NodeTrack, pltBeta=PlotBeta, SC=SpaceCharge), DemonNodePlacement)
        i += 1

    if(os.path.exists("BetaData_nosc.txt") and SpaceCharge == False):
        os.remove("betaData_nosc.txt")

    if(os.path.exists("BetaData_sc.txt") and SpaceCharge == True):
        os.remove("betaData_sc.txt")

    if(os.path.exists("BetaData.txt")):
        os.remove("betaData.txt")
        print "Cleared Beta Function Data"

#-----------------------------------------------------------------------------
# Initalize Lattice and Start Loop
#-----------------------------------------------------------------------------
print "Start Loop"

if RunLoops == False:
    nloops = 1

for i in range(nloops):
    # b.compress()
    print "Loop: "+str(i+1)

    if( (i == nloops-1 or nloops == 1) and DemonNodes == True):
        for Node in Lattice.getNodes():
            for ChildNode in Node.getChildNodes(DemonNodePlacement):
                if(ChildNode.getType() == "psDM"):
                    ChildNode.setLastLoop(True)
                    ChildNode.setDoBeta(True)

    Lattice.trackBunch(b, paramsDict)

    if PrintLoop == True:

        os.chdir(bf.DATpath)

        if(i < 10 and i < 100):
            index = "00" + str(i)
        elif(i < 100):
            index = "0" + str(i)
        else:
            index = str(i)

        fname = "loop" + index + ".dat"
        f = open(fname,'w+')

        for i in range(b.getSize()):
            coords = str(b.x(i))+" "+str(b.y(i))+" "+str(b.z(i))+'\n'
            # coords = coords+str(b.xp(i))+" "str(b.yp(i))+" "+str(b.zp(i))
            f.write(coords)
            # print(coords)
        f.close()
        os.chdir(bf.basepath)
print "End Loop"


b.dumpBunch("final.dat")


#-----------------------------------------------------------------------------
# File Management
#-----------------------------------------------------------------------------

#Move DemonNode Files
if DemonNodes == True and NodeTrack == True:
    for i in range(NodeCount):
        if(i < 10 and i < 100):
            index = "00" + str(i)
        elif(i < 100):
            index = "0" + str(i)
        else:
            index = str(i)

        name = "dm" + index + ".txt"

        bf.MoveFile(name, bf.TXTpath)



#Move data Files
files = ["initial.dat","final.dat"]
bf.MoveFile(files, bf.DATpath)

#-----------------------------------------------------------------------------
# Finished
#-----------------------------------------------------------------------------

print "Stop."
