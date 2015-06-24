# Simulation parameters

import numpy as np
import os

hallbar.outFile = "trans.npz"
hallbar.outDir = 'tmfsc_out/'

#hallbar.Bmin = 3#0.001
#hallbar.Bmax = 3
#hallbar.Vmin = 0.6#0.001
#hallbar.Vmax = 0.6
#hallbar.m = -1
#hallbar.NB = 2
#hallbar.NV = 2

hallbar.Bmin = 1.5
hallbar.Vmin = 0.3
hallbar.m = 1
hallbar.NB = 1
hallbar.NV = 1
hallbar.setupBias()
hallbar.printBias()

hallbar.drawGeom()
#hallbar.calcSingleTraject(shiftxy=(20.0,0.0), shiftth=-30.0/180.0*pi)
#hallbar.calcSingleTraject(contId = 3)
#hallbar.printTraj()
#hallbar.calcSingleTrans(saveTrajectory=True, contId=0)
#hallbar.calcSingleTrans(saveTrajectory=True, contId=3)
#hallbar.drawTrajectory()
#hallbar.animate()
#hallbar.showPlot()
hallbar.calcAllTrans()

# Save to a file
hallbar.save()




