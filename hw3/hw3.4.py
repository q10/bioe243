from hw3 import *
from math import exp, pi
from numpy import histogram
import pylab

## HW 3.4a STUFF

VV_sim = VelocityVerlet("LJ_108_1.txt")
num_timesteps = 1

time = []
MSD_table = []

#for t in range(1,num_timesteps+1):
#print "running step " + str(t) + ":"
t=5000
N = len(VV_sim.particles)
old_positions = []
    
for k in range(N):
    old_positions.append(VV_sim.particles[k].position())

for j in range(t):
    print "\trunning MD step " + str(j) + "..."
    VV_sim.MD_step()

mean_square_displacement = 0.0
for k in range(N):
    (x0 ,y0, z0) = old_positions[k]
    (x1, y1, z1) = VV_sim.particles[k].position()
    mean_square_displacement += ((x1-x0)**2)+((y1-y0)**2)+((z1-z0)**2)
mean_square_displacement /= (6*N*(t+1)*DT)

time.append((t+1)*DT)
MSD_table.append(mean_square_displacement)

pylab.clf()
pylab.plot(time, MSD_table)
#pylab.show()
#pylab.savefig('hw3.3a.y_histogram.png')



