# Simulation parameters

import numpy as np
import os
from math import sqrt, radians, degrees
from math import pi, tan, cos



# Geometry
# =============================================================================
def setupGeom(self):
    """Sets up the geometry."""
    factor = 1.0
    ly = 1.0*1000.0/factor          # width
    lx = 2*1000.0/factor          # length
    wc = 0.35*1000.0/factor         # contact width
    safeFactor = 0.0/factor   # safe factor for extending gate
    #xdc = 1.95*1000.0/factor  # half of horizontal distance between two contacts
    #ydc = ly / 2.0
    splitLen1 = 70                 # junction 1 Split nm
    splitLen2 = 70                  # junction 2 Split nm
    juncShift1 = 0
    juncShift2 = 0
    thr = 45                        # tilt angle in degree
    
    x1 = -lx/2;     y1 = -ly/2
    x2 = lx/2;      y2 = y1
    x3 = x2;        y3 = y1 + ly
    x4 = x1;        y4 = y3
    
    self.addPoint(x1, y1)
    self.addPoint(x2, y2)
    self.addPoint(x3, y3)
    self.addPoint(x4, y4)
    npts = self.addPoint(x1, y1) + 1
    self.mprint("\n-D- Number of vertices in geometry", npts)
    #add reflecting edges                                                  
    for ip in range(npts):                                                  
        self.addEdge(ip, (ip+1) % npts); 
    #add contacts
    for ip in range(1, npts):
        if ip == 3 or ip == 1 :
            #if ip == 2 or ip == 6 or ip == 12 or ip == 16:
            self.setEdgeType(ip, EDGE_ABSORB)
    thr = thr * pi/180                        # tilt angle in radian
    juncShift2 = juncShift2 * cos( thr )      # actual junction 2 shift due to tilt
    # Add the transmitting edge
    ipt1 = self.addPoint(x1 + lx/3.0 +juncShift1, -ly/2)
    ipt2 = self.addPoint(x1 + lx/3.0 +juncShift1, ly/2)

    ipt3 = self.addPoint(x1 + lx*2.0/3.0 +juncShift2 - ly*tan(thr)/2.0, -ly/2)
    ipt4 = self.addPoint(x1 + lx*2.0/3.0 +juncShift2 + ly*tan(thr)/2.0  ,  ly/2)
    self.addEdge(ipt1, ipt2, EDGE_TRANSMIT, splitLen1)
    self.addEdge(ipt3, ipt4, EDGE_TRANSMIT, splitLen2)
    # add gates
    self.addGate((x1-safeFactor,y1-safeFactor), (x1 + lx/3.0 +juncShift1, y1-safeFactor),
                 (x1 + lx/3.0 +juncShift1, y3+safeFactor), (x1-safeFactor, y3+safeFactor), 0.0)
    self.addGate((x1 + lx/3.0 +juncShift1, y1-safeFactor),
                 (x1 + lx*2.0/3.0 +juncShift2 - ly*tan(thr)/2.0, y1-safeFactor),
                 (x1 + lx*2.0/3.0 +juncShift2 + ly*tan(thr)/2.0, y3+safeFactor),
                 (x1 + lx/3.0 +juncShift1, y3+safeFactor), 0.0)
    self.addGate((x1 + lx*2.0/3.0 +juncShift2 - ly*tan(thr)/2.0, y1-safeFactor),
                 (x2+safeFactor, y2-safeFactor),
                 (x2+safeFactor, y3+safeFactor),
                 (x1 + lx*2.0/3.0 +juncShift2 + ly*tan(thr)/2.0, y3+safeFactor), -1.0)

HallBar.setupGeom = setupGeom

# Biasing scheme
# =============================================================================
def setupBias(self, Ef, B, V1, V2, V3, m = 1, singleResonance = True,
        Efmax = None, NEf=1, Bmax = None, NB = 1, V1max = None, NV1 = 1, 
                V2max=None, NV2=1, V3max = None, NV3=1):
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
    if NV3 == 1:
        VV3 = np.array([V3])
    else:
        assert V3max is not None, "V3max is None"
        VV3 = np.linspace(V3, V3max, NV3)

    self.bias.clear()
    self.bias.append(EEf, 'Ef')
    self.bias.append(BB, 'B')
    self.bias.append(VV1, 'V')
    self.bias.append(VV2, 'V')
    self.bias.append(VV3, 'V')

HallBar.setupBias = setupBias

hallbar.clear()
hallbar.outDir = 'APL_Paper/pnp_3Gate/Rgh_15/'
hallbar.cleanRun = True
hallbar.DebugLevel = 0
hallbar.verbosity = 1
hallbar.OccupationTol = 1E-4
hallbar.MinmNoInjection = 1500
hallbar.sim.TimeStep = 5.0/(1E6/1E-9)
#hallbar.EdgeRefRghnsEff = 1.0
#hallbar.sim.InjecAngleSpread = pi/3.0
EDGE_RGH_SPREAD = 15.0      # Degrees
JUNC_RGH_SPREAD = 10.0      # Degrees

hallbar.EdgRghAngleSpread = radians( EDGE_RGH_SPREAD )
hallbar.JuncRghAngleSpread = radians( JUNC_RGH_SPREAD )
hallbar.printParam()
hallbar.setupGeom()

EF = 0
n1 = 2.5E16
n2 = 2.5E16
V1 = 0.3
V2min = -V1
V2max = V1
Bmin = -0.25
Bmax = 0.25
hallbar.mprint("\n-D- V1:",V1)
hallbar.mprint("\n-D- V2min,max:",V2min,V2max)
hallbar.mprint("\n-D- Bmin,max:",Bmin,Bmax)
#hallbar.setupBias(B=0, V1=V1, V2=-V1, V3=V1, singleResonance=False)
hallbar.setupBias(Ef=-0.5, B=0, V1=V1, V2=V2min, V3=V1, V2max=V1, Efmax=0.5, 
                        NEf=3, NV2=3, singleResonance=False)

hallbar.printBiasList()

#hallbar.drawGeom(gateBorder = 1.0)
hallbar.enableCyclotronCalc()
#hallbar.enableDirectCalc()

#hallbar.calcSingleTraject(shiftth=-0.0/180.0*pi, contId = 0)
#hallbar.printTraj()

#hallbar.calcSingleTrans(saveTrajectory=True, contId=1)
hallbar.calcAllTrans(dl=5, nth=50, contId=1)

#hallbar.drawTrajectory(color='#343B9C', marker='', width=0.4)
#hallbar.animate()
#hallbar.showPlot()
#hallbar.saveTraj("struct2_nn_peak_theta_14p7_Sig_pi_5.png")
#hallbar.saveTraj("pnp_3Gate_Off_State.svg", dpi=600)

# Save to a file
hallbar.save()





