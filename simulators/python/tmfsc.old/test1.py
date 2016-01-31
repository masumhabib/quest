
import numpy as np
import matplotlib.pyplot as plt

# fermi velocity
nm = 1E-9   # nanometer
AA = 1E-10   # Angstrom

vf = 1E6
V = -0.15
EF = 0.0
Bz = 0.0
q = 1.60217657E10-19
e = -q

dt = 1E-15*1E-1
N = 10000

B = np.array([0.0, 0.0, Bz])
v = np.array([vf, 0.0, 0.0])

m = q*(EF - V)/vf**2.0

r = np.array([0,0,0])
ri = np.zeros((N, 3));
for i in range(N):
    F = e*np.cross(v, B)
    a = F/m
    #print (a)

    v = v + a*dt 
    r = r + v*dt
    ri[i, :] = r/AA


#print (ri)
fig = plt.figure()
axes = fig.add_subplot(111)
 
axes.plot(ri[:, 0], ri[:, 1], 'r-')
axes.set_aspect('equal', 'datalim')
plt.show()
