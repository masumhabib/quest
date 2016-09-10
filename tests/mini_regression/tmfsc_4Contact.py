# Simulation parameters

import numpy as np
import os
from math import sqrt
from math import pi
from math import exp as exponential

# Geometry
# =============================================================================
def setupGeom(self):
    """Sets up the geometry."""
    factor = 1
    ly = 3.9*1000.0/factor  # width
    lx = 7.9*1000.0/factor  # length
    wc = 0.3*1000.0/factor   # contact width
    wc2 = wc
    lc = 0.05*1000.0/factor   # contact length
    xdc = 1.95*1000.0/factor  # horizontal distance between two contacts
    ydc = ly / 2.0
    dg = 50/factor           # gate split
    xg = 5.0*1000.0/factor   # position of gate from left
    splitLen = 80.0 # gate Split nm
    juncShift = 0.0          # shift of physical junction
    x1 = -lx/2;     y1 = -ly/2
    x2 = -xdc;    y2 = y1;
    x3 = x2+2*xdc;  y3 = y2;
    x4 = lx/2;      y4 = y3;
    x5 = x4;        y5 = y4+ydc;
    x6 = x5;        y6 = y5+ydc;
    x7 = xdc;    y7 = y6;
    x8 = x7-2*xdc;  y8 = y7;
    x9 = -lx/2;    y9 = y8;
    x10 = x9;       y10 = y9-ydc;
    self.dc = xdc*2
    
    self.addPoint(x1, y1)
    self.addPoint(x2-wc2/2, y2)
    self.addPoint(x2-wc2/2, y2-lc)
    self.addPoint(x2+wc2/2, y2-lc)
    self.addPoint(x2+wc2/2, y2)
    self.addPoint(x3-wc/2, y3)
    self.addPoint(x3-wc/2, y3-lc)
    self.addPoint(x3+wc/2, y3-lc)
    self.addPoint(x3+wc/2, y3)
    self.addPoint(x4, y4)
    self.addPoint(x6, y6)
    self.addPoint(x7+wc/2, y7)
    self.addPoint(x7+wc/2, y7+lc)
    self.addPoint(x7-wc/2, y7+lc)
    self.addPoint(x7-wc/2, y7)
    self.addPoint(x8+wc/2, y8)
    self.addPoint(x8+wc/2, y8+lc)
    self.addPoint(x8-wc/2, y8+lc)
    self.addPoint(x8-wc/2, y8)
    self.addPoint(x9, y9)
    npts = self.addPoint(x1, y1) + 1
    self.mprint("\n-D- Number of vertices in geometry", npts)
    #add reflecting edges                                                  
    for ip in range(npts):                                                  
        self.addEdge(ip, (ip+1) % npts); 
    #add contacts
    for ip in range(1, npts):
        # Absorbing Side Contact
        if ip == 2 or ip == 6 or ip == 12 or ip == 16 or ip ==9 or ip == 19:
        # Without Absorbing Side Contact
        #if ip == 2 or ip == 6 or ip == 12 or ip == 16:
      		self.setEdgeType(ip, EDGE_ABSORB)
    # Add the transmitting edge
    ipt1 = self.addPoint(0+juncShift, -ly/2)
    ipt2 = self.addPoint(0+juncShift, ly/2)
    self.addEdge(ipt1, ipt2, EDGE_TRANSMIT, splitLen)
    # add gates
    self.addGate((-lx/2-lc,-ly/2-lc), (0+juncShift, -ly/2-lc),
        (0+juncShift, ly/2+lc), (-lx/2-lc, ly/2+lc), 0.0)
    self.addGate((0+juncShift,-ly/2-lc), (lx/2+lc, -ly/2-lc),
        (lx/2+lc, ly/2+lc), (0+juncShift, ly/2+lc), -1.0)
HallBar.setupGeom = setupGeom

# Biasing scheme
# =============================================================================
def setupBias(self, Ef, B, V1, V2, m = 1, singleResonance = True, 
        Efmax= None, NEf = 1, Bmax = None, NB = 1, V1max = None, NV1 = 1, V2max=None, NV2=1):
    """ Sets up the bias points """
    EEf = []
    BB = [] 
    VV1 = [] 
    VV2 = []
    if NEf == 1:   
        EEf = np.array([Ef])
    else:
        assert Efmax is not None, "Efmax is None"
        EEf = np.linspace(Ef, Efmax, NEf)  
 
    if NB == 1:   
        BB = np.array([B])
    else:
        assert Bmax is not None, "Vmax is None"
        BB = np.linspace(B, Bmax, NB)
    if singleResonance == True:
        B0 = abs(self.EF-V1)/vf/nm/self.dc*2.0
        BB = m*np.array([B0])

    if NV1 == 1:
        VV1 = np.array([V1])
    else:
        assert V1max is not None, "V1max is None"
        VV1 = np.linspace(V1, V1max, NV1)
    if NV2 == 1:
        VV2 = np.array([V2])
    else:
        assert V2max is not None, "V2max is None"
        VV2 = np.linspace(V2, V2max, NV2)

    self.bias.clear()
    self.bias.append(EEf, 'Ef')
    self.bias.append(BB, 'B')
    self.bias.append(VV1, 'V')
    self.bias.append(VV2, 'V')

HallBar.setupBias = setupBias

hallbar.clear()
## Simulation Parameters
# output directory
hallbar.outDir = 'dummyOut/'
hallbar.DebugLevel = 0
hallbar.verbosity = 1
hallbar.OccupationTol = 1E-4
# Minimum number of electron injection
hallbar.MinmNoInjection = 100
# If EdgeRefRghnsEff is set to any other value than 1.0,
# it will activate crude model
hallbar.EdgeRefRghnsEff = 1.0
# If InjecAngleSpread is not set, cosine distribution will be used
hallbar.InjecAngleSpread = pi/15
hallbar.setupGeom()

hallbar.mprint("\n-D- OccupationTol", hallbar.OccupationTol)
hallbar.mprint("\n-D- Injection spread angle:", hallbar.sim.InjecAngleSpread)

## Biasing Parameters
V1 = 0.0928
V2min = 0.0928
V2max = 0.0928
Bmin = 0.0476
Bmax = 0.0476
hallbar.mprint("\n-D- V1:", V1)
hallbar.mprint("\n-D- V2min,max:", V2min, V2max)
hallbar.mprint("\n-D- Bmin,max:", Bmin, Bmax)

hallbar.setupBias(Ef=0.0, B=Bmin, V1=V1, V2=V2min, Efmax= 2.0, NEf = 10, Bmax=Bmax, V2max=V2max, NB=1, NV2=1, singleResonance=False)
# for setting custom gate color
# gateColors = np.array( [[1, 1, 1], [0.816, 0.812, 0.792]] )
# hallbar.drawGeom(gateBorder = 1.0, gateColors = gateColors)
hallbar.printBiasList()

hallbar.drawGeom(gateBorder = 1.0)
hallbar.enableCyclotronCalc() 
# Single Bias point calculation - for seeing trajectory
hallbar.calcSingleTrans(saveTrajectory=True, contId=0)
hallbar.drawTrajectory(color='#343B9C', marker='', width=0.3)
# All Bias point calculation - no trajectory
#hallbar.calcAllTrans(dl=10, nth=20, contId=%CONTID%)
# saving trajectory image
hallbar.saveTraj("trajectory.png", dpi=600)

# Save to a file all bias point or single bias point data
#hallbar.save()





