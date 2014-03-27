"""  
 Author: K M Masum Habib <masum.habib@virginia.edu>
 Last update: 03/26/2014
"""

import numpy as np

"""
    Returns numpy.ndarray of regularly spaced grid. 
    It is similar to numpy.linspace except that this function
    takes minimum, maximum and delta as the arguments.
"""
def linspace(min, max, d):
    n = int(abs(max - min)/d + 1)
    return np.linspace(min, max, n)
