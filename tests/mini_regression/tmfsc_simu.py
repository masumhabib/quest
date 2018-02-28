# Simulation parameters

import numpy as np
import os

hallbar.outFile = "trans.npz"
hallbar.outDir = 'tmfsc_out/'

#hallbar.setupBias(0.0, 0.001, m=1, singleResonance=False,
#        Bmax=0.0, NB=2, Vmax=0.001, NV=2)

hallbar.setupBias(Ef = 0.0, B=1.5, V=0.15)
hallbar.printBiasList()
#hallbar.sim.TimeStep = 1/(1E6/1E-9)

# FOR GATES
lx = 500.0		# length
ly = 500.0		# width
clx = 50.0		# contact width
cly = 1.0		# contact length
coffx = lx/10
# single gate
#hallbar.addGate((-lx/2,-ly/2-cly*50), (lx/2, -ly/2-cly*50),
#                (lx/2, ly/2), (-lx/2, ly/2))
# double gate, pn junction
#hallbar.addGate((-lx/2,-ly/2-cly*50), (0, -ly/2-cly*50),
#                (0, ly/2), (-lx/2, ly/2), 1)
#hallbar.addGate((0,-ly/2-cly*50), (lx/2, -ly/2-cly*50),
#                (lx/2, ly/2), (0, ly/2), -1)

#hallbar.enableDirectCalc()
hallbar.enableCyclotronCalc()
hallbar.drawGeom()
#hallbar.calcSingleTraject(shiftxy=(20.0,0.0), shiftth=-30.0/180.0*pi)
#hallbar.calcSingleTraject(shiftth=-30.0/180.0*pi)
#hallbar.calcSingleTraject(shiftxy=(30.0,200.0))
hallbar.calcSingleTraject(contId = 0)
#hallbar.printTraj()
#hallbar.calcSingleTrans(saveTrajectory=True, contId=0)
#hallbar.calcSingleTrans(saveTrajectory=True, contId=0)
hallbar.drawTrajectory()
#hallbar.animate()
hallbar.showPlot()
#hallbar.calcAllTrans()

# Save to a file
#hallbar.save()




