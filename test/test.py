#!/usr/bin/python

import qmicad as qm
import numpy as np

x = np.array([[1., 2., 3.],[4., 5., 6]])
qm.test_npy2mat(x)
qm.test_npy_to_mat(x)

y = qm.test_mat_to_npy()
qm.test_npy2mat(y)
qm.test_npy_to_mat(y)

