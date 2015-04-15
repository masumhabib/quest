#!/usr/bin/python

import qmicad as qm
import numpy as np 

A = np.zeros( (10, 1), dtype=np.double)
B = np.zeros( (5, 1), dtype=np.double)
temp = qm.potential.poissonPot(A, B)
print temp