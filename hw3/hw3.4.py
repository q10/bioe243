from hw3 import *
from math import exp, pi
import sys

#from numpy import histogram
#import pylab

## HW 3.4a STUFF


if len(sys.argv) != 2:  # the program name and the two arguments
  # stop the program and print an error message
  sys.exit("Must provide a number [0-3]")

# select pressure and temperature pairs
RHO = [0.8442, 0.8442, 0.88, 0.92]
TARGET_TEMP = [0.742, 0.9, 0.8, 0.742]

# change these variables
SELECT = int(sys.argv[1])

# initialize system
VV_sim = VelocityVerlet("LJ_108_1.txt", density=RHO[SELECT], temp=TARGET_TEMP[SELECT])
timesteps=30000
N = len(VV_sim.particles)
old_positions = []
avg_VV_T = 0.0

# record positions at t=0
for k in range(N):
    old_positions.append(VV_sim.particles[k].position())

# run simulation
for j in range(timesteps):
    #print "\trunning MD step " + str(j) + "..."
    VV_sim.MD_step()
    avg_VV_T += VV_sim.system_kinetic_energy

# calculate mean square displacement and average temperature
mean_square_displacement = 0.0
for k in range(N):
    (x0 ,y0, z0) = old_positions[k]
    (x1, y1, z1) = VV_sim.particles[k].position()
    mean_square_displacement += ((x1-x0)**2)+((y1-y0)**2)+((z1-z0)**2)

mean_square_displacement /= (6*N*timesteps*DT)
avg_VV_T = (2*(avg_VV_T/timesteps))/((3*N)-3)

# write results to files
if SELECT == 0:
    f = open('hw3.4.table_row_1.results','w')
elif SELECT == 1:
    f = open('hw3.4.table_row_2.results','w')
elif SELECT == 2:
    f = open('hw3.4.table_row_3.results','w')
else:
    f = open('hw3.4.table_row_4.results','w')

print >>f, "Input rho* = " + str(RHO[SELECT]) + "\nInput temp* = " + str(TARGET_TEMP[SELECT]) + "\nAverage temp* = " + str(avg_VV_T) + "\nSimulated D* = " + mean_square_displacement
f.close()
