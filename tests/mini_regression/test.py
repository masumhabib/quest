#!/usr/bin/python

import quest as qm
import numpy as np

x = np.array([[1., 2., 3.],[4., 5., 6]])
qm.test_npy_to_mat(x)
#qm.test_npy2mat(x)


y = qm.test_mat_to_npy()
qm.test_npy_to_mat(y)
#qm.test_npy2mat(y)

yy = qm.test_mat_to_npy()
qm.test_npy_to_mat(yy)
#qm.test_npy2mat(yy)


#xx = np.array([[1.+1.j, 2.+2.j, 3.+3.j],[4.+4.j, 5.+5.j, 6.+6.j]])
#qm.test_npy_to_mat(xx)
#qm.test_npy2mat(xx)
