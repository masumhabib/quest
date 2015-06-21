# Simulation parameters

import numpy as np
import os

hallbar.out_file = "trans.npz"
hallbar.out_dir = 'tmfsc_out/'

hallbar.num_contacts = 4

hallbar.Bmin = 3#0.001
hallbar.Bmax = 3
hallbar.Vmin = 0.6#0.001
hallbar.Vmax = 0.6
hallbar.m = -1
hallbar.NB = 2
hallbar.NV = 2

hallbar.Bmin = 1.5
hallbar.Vmin = 0.3
hallbar.m = 1
hallbar.NB = 1
hallbar.NV = 1

hallbar.setup_geom()
hallbar.setup_bias()
hallbar.print_bias()
hallbar.draw_geom()
#hallbar.calc_single_traject(shiftxy=(20.0,0.0), shiftth=-30.0/180.0*pi)
#hallbar.print_traject()
hallbar.calc_single_trans(saveTrajectory=True)
hallbar.draw_trajectory()
#hallbar.animate()
hallbar.show_plot()
#hallbar.calc_all_traject()
#hallbar.calc_all_trans()

# Save to a file
#hallbar.save()




