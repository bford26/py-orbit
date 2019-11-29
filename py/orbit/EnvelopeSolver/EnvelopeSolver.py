
import envelope_solver as ES

# Imports
import sys
import os
import math
import numpy as np


from orbit.teapot.teapot import NodeTEAPOT, TiltTEAPOT

from orbit.teapot_base import TPB
# import the function that creates multidimensional arrays
from orbit.utils import orbitFinalize
# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker
# import the MAD parser to construct lattices of TEAPOT elements.
from orbit.parsers.mad_parser import MAD_Parser, MAD_LattElement
# import the MADX parser to construct lattices of TEAPOT elements.
from orbit.parsers.madx_parser import MADX_Parser, MADX_LattElement
# import aperture
from aperture import Aperture
# monitor
from bunch import BunchTwissAnalysis


class BasicSolverNode(NodeTEAPOT):

    def __init__(self, dL=0.0, Disp=0.0 ):
        AccNodeBunchTracker.__init__(self,name="no name BS Node")
        self.setType("Basic_EnvSolver")
        self.dLen = dL
        self.Dispersion = Disp
        self.SC = 0.0

    def setDispersion(self, Disp=0.0):
        self.Dispersion = Disp

    def setDLength(self, len=0.0):
        self.dLen = len

    def setSpaceCharge(self, sc=0.0):
        self.SC = sc

    def track(self, paramsDict):
        EnvBunch = paramsDict["bunch"]
        EmitX = paramsDict["EmitX"]
        EmitY = paramsDict["EmitY"]

        if paramsDict.has_key("spacecharge"):
            SpaceCharge = paramsDict["spacecharge"]
        else:
            SpaceCharge = self.SC

        if paramsDict.has_key("dlength"):
            dLen = paramsDict["dlength"]
        else:
            dLen = self.dLen

        if paramsDict.has_key("dispersion"):
            Dispersion = paramsDict["dispersion"]
        else:
            Dispersion = self.Dispersion

        ES.BasicSolver(EnvBunch, EmitX, EmitY, SpaceCharge, Dispersion, dLen)

class TiltSolverNode(NodeTEAPOT):
    """
    Envelope Solver Node for tilting beam envelope
    """
    def __init__(self, spacecharge=0.0, length=0.0):
        AccNodeBunchTracker.__init__(self, name="no name TiltSol Node")
        self.setType("Tilt_EnvSolver")
        self.dLen = length
        self.SC = spacecharge

    def track(self, paramsDict):
        EnvBunch = paramsDict["bunch"]

        if paramsDict.has_key("spacecharge"):
            SpaceCharge = paramsDict["spacecharge"]
        else:
            SpaceCharge = self.SC

        if paramsDict.has_key("dlength"):
            dLen = paramsDict["dlength"]
        else:
            dLen = self.dLen

        ES.LinacTiltSolver(EnvBunch, SpaceCharge, dLen)











#eof
