# Simulation parameters

import numpy as np
import os

hallbar.outFile = "trans.npz"
hallbar.outDir = 'tmfsc_out/'

hallbar.Bmin = 0.0#0.001
hallbar.Bmax = 0.0
hallbar.Vmin = 0.001#0.001
hallbar.Vmax = 0.001
hallbar.m = 1
hallbar.NB = 2
hallbar.NV = 2

#hallbar.Bmin = 1.5
#hallbar.Vmin = 0.15
#hallbar.m = 1
#hallbar.NB = 1
#hallbar.NV = 1
hallbar.setupBias()
hallbar.printBias()
hallbar.sim.TimeStep = 1/(1E6/1E-9)

hallbar.enableDirectCalc()
#hallbar.enableCyclotronCalc()
hallbar.drawGeom()
#hallbar.calcSingleTraject(shiftxy=(20.0,0.0), shiftth=-30.0/180.0*pi)
#hallbar.calcSingleTraject(shiftxy=(30.0,200.0))
#hallbar.calcSingleTraject(contId = 0)
#hallbar.printTraj()
#hallbar.calcSingleTrans(saveTrajectory=True, contId=0)
hallbar.calcSingleTrans(saveTrajectory=True, contId=0)
hallbar.drawTrajectory()
#hallbar.animate()
hallbar.showPlot()
#hallbar.calcAllTrans()

# Save to a file
#hallbar.save()




