from hw3 import *
from common import *
from math import exp, pi
from numpy import histogram
import pylab
'''
## HW3.3a STUFF

r = BoxMueller(50232364960420517)
x_histogram, y_histogram = [], []
for i in range(0,10000):
    (x,y) = r.generate()
    x_histogram.append(x)
    y_histogram.append(y)

x_histogram = histogram(x_histogram, bins=70)
y_histogram = histogram(y_histogram, bins=70)

# save data
f = open('hw3.3a.x_histogram.dat','w')
print >>f, x_histogram
f.close()
f = open('hw3.3a.y_histogram.dat','w')
print >>f, y_histogram
f.close()

# plot data and save
pylab.clf()
(n, bins) = x_histogram
pylab.plot(.5*(bins[1:]+bins[:-1]), n)
pylab.savefig('hw3.3a.x_histogram.png')

pylab.clf()
(n, bins) = y_histogram
pylab.plot(.5*(bins[1:]+bins[:-1]), n)
pylab.savefig('hw3.3a.y_histogram.png')
'''


## HW3.3c STUFF

steps = 5000

VV_sim = VelocityVerlet("LJ_108_1.txt")
FE_sim = ForwardEuler("LJ_108_1.txt")

avg_VV_KE, avg_FE_KE = 0.0, 0.0
time = [0.0]
VV_KE = [VV_sim.system_kinetic_energy]
VV_PE = [VV_sim.system_potential_energy]
VV_TE = [VV_sim.system_total_energy]
FE_KE = [FE_sim.system_kinetic_energy]
FE_PE = [FE_sim.system_potential_energy]
FE_TE = [FE_sim.system_total_energy]

for i in range(steps):
    print "step " + str(i) + " running..."
    VV_sim.MD_step()
    FE_sim.MD_step()

    avg_VV_KE += VV_sim.system_kinetic_energy
    avg_FE_KE += FE_sim.system_kinetic_energy

    time.append(VV_sim.time)
    VV_KE.append(VV_sim.system_kinetic_energy)
    VV_PE.append(VV_sim.system_potential_energy)
    VV_TE.append(VV_sim.system_total_energy)
    FE_KE.append(FE_sim.system_kinetic_energy)
    FE_PE.append(FE_sim.system_potential_energy)
    FE_TE.append(FE_sim.system_total_energy)

avg_VV_KE /= steps
avg_FE_KE /= steps
avg_VV_T = (2*avg_VV_KE)/((3*len(VV_sim.particles))-3)
avg_FE_T = (2*avg_FE_KE)/((3*len(FE_sim.particles))-3)


f = open('hw3.3c.results','w')
print >>f, "Average temperature in Velocity Verlet: " + str(avg_VV_T)
print >>f, "Average temperature in Forward Euler: " + str(avg_FE_T) + "\n\n"
for i in range(steps+1):
    print >>f, str(time[i]) + "\t" + str(VV_KE[i]) + "\t" + str(VV_PE[i]) + "\t" + str(VV_TE[i]) + "\t" + str(FE_KE[i]) + "\t" + str(FE_PE[i]) + "\t" + str(FE_TE[i]) + "\n"
f.close()

pylab.clf()
pylab.plot(time, VV_KE, color='blue')
pylab.plot(time, VV_PE, color='red')
pylab.plot(time, VV_TE, color='green')
pylab.savefig('hw3.3c.VV.png')

pylab.clf()
pylab.plot(time, FE_KE, color='blue')
pylab.plot(time, FE_PE, color='red')
pylab.plot(time, FE_TE, color='green')
pylab.savefig('hw3.3c.FE.png')
