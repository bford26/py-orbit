from mpl_toolkits.mplot3d import Axes3D

import numpy as np
import matplotlib.pyplot as plt

import math
import sys
import os, glob
import shutil

class FileFunctions:

    def __init__(self):
        self.basepath = None
        self.TXTpath = None
        self.DATpath = None
        self.PLOTpath = None

        self.xlim = 0.1
        self.ylim = 0.1

    def setTextPath(self, TextPath):
        self.TXTpath = TextPath
    def setDatPath(self, DatPath):
        self.DATpath = DatPath
    def setPlotPath(self, PlotPath):
        self.PLOTpath = PlotPath
    def setBasePath(self, BasePath):
        self.basepath = BasePath

    def MoveFile(self, File, Destination, BasePath=None):
        """
        This function will move files to prefered directories
        """
        if BasePath == None:
            if self.basepath == None:
                print "Please Set BasePath For File Directory"
            else:
                BasePath = self.basepath

        if type(File) == type("str"):
            name = File
            if(os.path.exists(Destination+name) == False):
                shutil.move(BasePath+name,Destination+name)
            else:
                os.remove(Destination+name)
                shutil.move(BasePath+name,Destination+name)
        else:
            for name in File:
                if(os.path.exists(Destination+name) == False):
                    shutil.move(BasePath+name,Destination+name)
                else:
                    os.remove(Destination+name)
                    shutil.move(BasePath+name,Destination+name)

    def FindFiles(self, ext=None):
        """
        With Extension set to None the returned list will be .txt files
        """
        if ext == None:
            ext = ".txt"

        list = []

        if ext == ".txt":
            path = self.TXTpath
        elif ext == ".dat":
            path = self.DATpath
        else:
            print "Extension not recognized"

        os.chdir(path)
        for file in glob.glob('*'+ext):
            if file.endswith(ext):
                list.append(file)

        os.chdir(self.basepath)
        list.sort()
        return list

    def ReadData(self, FileName, Delim=None, LineLimit=None):
        """
        Read Data from files giving in particular format:
        ' ' == deliminator , but has option for more
        % commenting
        x x2 y y2 z z2 style format
        """

        if FileName.endswith(".txt"):
            fpath = self.TXTpath
        elif FileName.endswith(".dat"):
            fpath = self.DATpath

        if Delim == None:
            Delim = ' '
        f = open(fpath + FileName, 'r')
        xdata = []
        ydata = []
        zdata = []
        xpd = []
        ypd = []
        zpd = []
        for i,str in enumerate(f.readlines()):
            str0 = str
            if str0[0] != '%':
                p = str0.split(Delim)
                #print(p)
                if FileName[0]=='l':
                    #0 1 2
                    #x y z
                    x = p[0]
                    y = p[1]
                    z = p[2]
                else:
                    #0 1  2 3  4 5
                    #x xp y yp z zp
                    x = p[0]
                    y = p[2]
                    z = p[4]
                    xp = p[1]
                    yp = p[3]
                    zp = p[5]
                    xpd.append(float(xp))
                    ypd.append(float(yp))
                    zpd.append(float(zp))
                xdata.append(float(x))
                ydata.append(float(y))
                zdata.append(float(z))

                if(type(LineLimit) == type(2) and LineLimit == i-1):
                    break

        f.close()
        xdata = np.array(xdata)
        ydata = np.array(ydata)
        zdata = np.array(zdata)
        xpd = np.array(xpd)
        ypd = np.array(ypd)
        zpd = np.array(zpd)
        return [xdata, xpd, ydata, ypd, zdata, zpd]

    def clearDM(self):
        files = self.FindFiles()
        # if len(files) >
        for file in files:
            os.remove(self.TXTpath + file)
        print "DM files removed"

    def clearLoops(self):
        files = self.FindFiles(".dat")
        # files.remove("initial.dat")
        # files.remove("final.dat")
        for file in files:
            os.remove(self.DATpath + file)
        print "Loop files removed"

    def clearPlots(self):

        path = self.PLOTpath
        SP = "SpecificPlots/"
        plotpaths = [ "", "dmplots/","loops/",SP+"xxp/" , SP+"xyp/" , SP+"yxp/" , SP+"yyp/"]

        for i in range(len(plotpaths)):
            os.chdir(path + plotpaths[i])
            for file in glob.glob('*.png'):
                os.remove(path+plotpaths[i]+file)
        print "Plots removed"

    def RemoveAllData(self):
        self.clearDM()
        self.clearLoops()

class BeamFunctions(FileFunctions):
    """
    Functions for beam physics and values
    """
    def setPlotLimits(self, xlim, ylim):
        self.xlim = xlim
        self.ylim = ylim

    def setSpaceCharge(self, SpaceCharge=None):
        if SpaceCharge == None:
            self.spacecharge = False
        else:
            self.spacecharge = SpaceCharge

    def eps_rms(self, u, up):
        """
        Statistical Emittance for u, either x or y coordinate
        """
        sigma_u = np.std(u)
        sigma_up = np.std(up)
        sigma_u_up = self.sig2(u,up)

        sigma_u_u = sigma_u*sigma_u
        sigma_up_up = sigma_up*sigma_up
        eps = math.sqrt((sigma_u_u * sigma_up_up) - (sigma_u_up*sigma_u_up))
        return eps

    def beta_stat(self, u, up):
        return np.std(u)*np.std(u)/(self.eps_rms(u,up))

    def gamma_stat(self,u, up):
        return np.std(up)*np.std(up)/(self.eps_rms(u,up))

    def alpha_stat(self,u, up):
        return math.sqrt(self.beta_stat(u,up)*self.gamma_stat(u,up) - 1)

    def twiss_stat(self,u, up):
        alpha = self.alpha_stat(u,up)
        beta = self.beta_stat(u,up)
        gamma = self.gamma_stat(u,up)
        return (alpha, beta, gamma)

    def sig2(self,u,up):
        sum1 = 0
        N = len(u)
        uavg = np.average(u)
        upavg = np.average(up)
        for i in range(N):
            sum1 = sum1 + (u[i] - uavg) * (up[i] - upavg)
        return abs(sum1 / N)

    def PlotBetaFunction(self, TrackingFilesDir=None, ext=None):
        """
        Plots the Beta function of of the beam. Defaults ext to .txt
        """
        if ext == None:
            ext = '.txt'

        if TrackingFilesDir == None and self.TXTpath == None:
            print "No Directory to get data from..."
            return

        if self.spacecharge == True:
            sc = " - SpaceCharge"
            sc1 = "SC.png"
        else:
            sc = " - No SpaceCharge"
            sc1 = "NO_SC.png"

        FigureName = "BetaPlot_" + sc1

        FileList = self.FindFiles(ext)

        StatBetaX = np.zeros((len(FileList),1))
        StatBetaY = np.zeros((len(FileList),1))

        i = 0
        for name in FileList:
            d = self.ReadData(self.TXTpath + name)
            StatBetaX[i] = self.beta_stat(d[0],d[1])
            StatBetaY[i] = self.beta_stat(d[2],d[3])
            i += 1

        plt.subplot(2, 1, 1)
        plt.plot(StatBetaX)
        plt.ylabel("Amp BetaX")
        plt.title("Beta Function"+sc)

        plt.subplot(2, 1, 2)
        plt.plot(StatBetaY)
        plt.xlabel("Nodes")
        plt.ylabel("Amp BetaY")

        plt.savefig(FigureName)
        plt.close()

    def plotANY(self, file, phasespace=None):

        """
        This function will take the phase space combination string
        and make the plot for the given file
        """
        if len(file) == 0:
            print "No Files Found"
            return

        if file.endswith(".txt"):
            path = self.TXTpath
            plt_path = self.PLOTpath+"dmplots/"
        elif file.endswith(".dat"):
            path = self.DATpath
            plt_path = self.PLOTpath+"loops/"


        if self.spacecharge ==True:
            sc = "SpaceCharge"
        else:
            sc = "No SpaceCharge"

        data = self.ReadData(file)

        if phasespace == None:
            phasespace = "xy"
        if phasespace == "xy":
            xaxis = data[0]
            yaxis = data[2]
            xla = "x [m]"
            yla = "y [m]"
        elif phasespace == "xyp":
            xaxis = data[0]
            yaxis = data[3]
            xla = "x [m]"
            yla = "yp [mrad]"
            plt_path = self.PLOTpath + "SpecificPlots/" + phasespace + '/'
        elif phasespace == "yxp":
            xaxis = data[2]
            yaxis = data[1]
            xla = "y [m]"
            yla = "xp [mrad]"
            plt_path = self.PLOTpath + "SpecificPlots/" + phasespace + '/'
        elif phasespace == "xxp":
            xaxis = data[0]
            yaxis = data[1]
            xla = "x [m]"
            yla = "xp [mrad]"
            plt_path = self.PLOTpath + "SpecificPlots/" + phasespace + '/'
        elif phasespace == "yyp":
            xaxis = data[2]
            yaxis = data[3]
            xla = "y [m]"
            yla = "yp [mrad]"
            plt_path = self.PLOTpath + "SpecificPlots/" + phasespace + '/'

        if(file[0] == 'i' or file == "dm000xy.txt"):
            self.xlim = max(xaxis)
            self.ylim = max(yaxis)

        index = file.rfind('.')
        file_type = file[index:]
        name = file[0:index]

        FigName = name+phasespace+".png"

        os.chdir(self.basepath)

        plt.figure(dpi = 200)
        plt.scatter(xaxis,yaxis)
        plt.title(name+phasespace+" "+sc)

        plt.xlim(-1*self.xlim,self.xlim)
        plt.ylim(-1*self.ylim,self.ylim)

        plt.xlabel(xla)
        plt.ylabel(yla)
        plt.savefig(FigName)
        plt.close()

        self.MoveFile(FigName,plt_path)

    def FastBetaPlot(self):

        os.chdir(self.basepath)

        f = open("BetaData.txt")
        if f.closed:
            print 'file is closed'
        # else:
        #     print 'file opened'

        nElements = sum(1 for line in open('BetaData.txt'))

        betaX = np.zeros((nElements,1))
        betaY = np.zeros((nElements,1))

        i=0
        for line in f.readlines():
            betas = line.split()
            betaX[i] = float(betas[0])
            betaY[i] = float(betas[1])
            i+=1

        f.close()

        if self.spacecharge == True:
            sc = " - SpaceCharge"
            sc1 = "SC.png"
        else:
            sc = " - No SpaceCharge"
            sc1 = "NO_SC.png"

        FigureName = "BetaPlot_" + sc1


        plt.figure()

        plt.plot(betaX)
        plt.plot(betaY)

        plt.xlabel("Nodes")
        plt.ylabel("Amp Beta")
        plt.title("Beta Function"+sc)
        plt.legend(["Beta X","Beta Y"])

        plt.savefig(FigureName)
        plt.close()

    def PosterPlot(self, nParts, LatticeLength, LatType):

        FileList = self.FindFiles(".txt")

        if type(FileList) == type("string"):
            print "Found Only One File"
            return
        else:
            nElems = len(FileList)

        xd = np.zeros(shape=(nParts,nElems))
        yd = np.zeros(shape=(nParts,nElems))
        zd = np.zeros(shape=(nParts,nElems))

        for j,filename in enumerate(FileList):
            data = self.ReadData(filename, Delim=None,LineLimit=nParts)
            xd[:,j] = data[0][:,0]
            yd[:,j] = data[2][:,0]
            zd[:,j] = data[4][:,0]

        if LatType == "LINAC":
            for i,j in zip(range(nElements), range(nParts)):
                Z = np.linspace(0,LatticeLength,nElems)
                for i in range(nParts):
                    for j in range(nElems):
                        zd[i,j] = zd[i,j] + Z[j]
        elif LatType == "RING":
            r = LatticeLength / (2 * math.pi)
            Theta = np.linspace(0,2 * math.pi,nElems)
            for i in range(nParts):
                for j in range(nElems):
            #         zd[i,j] = zd[i,j] + Z[j]
                    zd[i,j] = zd[i,j] + r * math.sin(Theta[j])
                    xd[i,j] = xd[i,j] - r * math.cos(Theta[j])

class EnvelopeFunctions(FileFunctions):

    def PlotAmplitudes(self):

        FileList = self.FindFiles(".txt")

        a = np.zeros((len(FileList),1))
        b = np.zeros((len(FileList),1))
        # ap = np.zeros((len(FileList),1))
        # bp = np.zeros((len(FileList),1))
        # Sig = np.zeros((len(FileList),1))

        for i,file in enumerate(FileList):
            f = open(self.TXTpath + file)
            for line in f.readlines():
                if line[0] != '%':
                    data = line.split()
                    a[i] = float(data[0])
                    b[i] = float(data[2])
                    # ap[i] = float(data[1])
                    # bp[i] = float(data[3])


        plt.subplot(2,1,1)
        plt.plot(a)
        plt.ylabel("A Amp")
        # plt.xlim(-self.xlim,self.xlim)

        plt.subplot(2,1,2)
        plt.plot(b)
        plt.ylabel("B Amp")
        plt.xlabel("Node")

        plt.savefig("AmpPlot.png")
        plt.close()

    def setLimits(self,xlim,ylim):
        self.xlim = xlim
        self.ylim = ylim


























#eof
