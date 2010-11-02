from math import exp, sqrt
from common import Ran3, BoxMueller, Particle, System
from hw41 import LocalMin_LJ

# default parameters
SEED = 3231582028467
SEED_2 = 845340327839343
SEED_3 = 324823012541
INIT_TEMP = 1.5
RHO = 0.01
TOLERANCE = 0.001
DT = 0.001

# define optimization
class SimulatedAnneal_LJ(System):
    def __init__(self, filename, temp=INIT_TEMP, density=RHO, time=DT, r_seed=SEED, schedule="geometric"):
        System.__init__(self, filename, temp, density)
        
        # take locally minimized solution from hw4.1 and copy over to start state
        temp_LocalMin_LJ = LocalMin_LJ(filename)
        temp_LocalMin_LJ.run()
        temp_positions_table = temp_LocalMin_LJ.collect_positions()
        self._restore_old_positions(temp_positions_table)
        self._update_forces()
        self._update_energies()
        
        self._initialize_velocities()
        self.time, self.step_time, self.step_time_2 = 0, time, time**2.0
        self.rand = Ran3(r_seed)
        self.cooling_schedule = schedule

    # set initial velocities according to Gaussian distribution
    def _initialize_velocities(self):
        vrand_1 = BoxMueller(SEED_2)
        vrand_2 = BoxMueller(SEED_3)

        for particle in self.particles:
            tmp = vrand_1.generate()
            particle.set_velocity(tmp[0], tmp[1], vrand_2.generate()[0])

    def _rescale_velocities(self):
        # T* = Tk_b / e ==> T* = sum(mv^2)/3N
        current_temp = 0.0
        for particle in self.particles:
            current_temp += particle.m * ((particle.vx**2) + (particle.vy**2) + (particle.vz**2))
        current_temp /= 3*self.N
        
        if abs(current_temp-self.temperature) > 0.001:
            gamma = sqrt(self.temperature/current_temp)
            for particle in self.particles:
                particle.vx, particle.vy, particle.vz = particle.vx*gamma, particle.vy*gamma, particle.vz*gamma

    # MD step runs on Velocity Verlet Algorithm
    def MD_step(self):
        #self.time += DT
        for first in self.particles:           
            # second terms of Taylor expansion, for less computing in next 2 steps
            sx = first.fx/(2.0*first.m)
            sy = first.fy/(2.0*first.m)
            sz = first.fz/(2.0*first.m)
            
            # update position         
            dx = (first.vx*self.step_time) + (sx*self.step_time_2)
            dy = (first.vy*self.step_time) + (sy*self.step_time_2)
            dz = (first.vz*self.step_time) + (sz*self.step_time_2)
            first.add_position(dx, dy, dz)

            # update velocity at half step with current forces
            dvx = sx*self.step_time
            dvy = sy*self.step_time
            dvz = sz*self.step_time
            first.add_velocity(dvx, dvy, dvz)

            # reset force
            first.fx, first.fy, first.fz = 0.0, 0.0, 0.0
            
        # update forces and energies
        self._update_forces()
        self._update_energies()
        
        # update velocity at second half step with new forces
        for first in self.particles:
            dvx = (first.fx*self.step_time)/(2*first.m)
            dvy = (first.fy*self.step_time)/(2*first.m)
            dvz = (first.fz*self.step_time)/(2*first.m)
            first.add_velocity(dvx, dvy, dvz)

        # rescale velocities to fit target temperature
        self._rescale_velocities()

    def save_old_system(self):
        return(self.collect_positions(), self.collect_velocities(), self.collect_forces(), self.collect_energies())

    def _restore_old_system(self, old_system_table):
        self._restore_old_positions(old_system_table[0])
        self._restore_old_velocities(old_system_table[1])
        self._restore_old_forces(old_system_table[2])
        self._restore_old_energies(old_system_table[3])
    
    def run(self):
        while self.temperature > 0.001:
            # run cooling schedule
            if self.cooling_schedule == "linear":
                self.temperature -= 0.01
            # else the cooling schedule is geometric
            else:
                self.temperature *= 0.99
                
            old_system_table = self.save_old_system()

            # Run MD moves
            steps = int(10 / self.temperature)
            for k in range(steps):
                self.MD_step()
            
            # averages
            energy_change = self.system_potential_energy() - sum(old_system_table[3])/2.0
            if energy_change > 0.0 and self.rand.generate() >= exp(-energy_change/self.temperature):
                # Revert back to state prior to simulating at current temperature,
                # since the new state has a higher energy AND rand is NOT < exp(-dE)
                self._restore_old_system(old_system_table)
            # print self.system_potential_energy()
                

###########################
# RUN CODE

f = open('hw4.2.results','w')
seed = Ran3(1249834071)

for i in range(2):
    filename = 'coord.' + str(i+1)
    for j in range(3):
        a = SimulatedAnneal_LJ(filename, r_seed=int(1E12*seed.generate()), schedule="geometric")
        b = SimulatedAnneal_LJ(filename, r_seed=int(1E12*seed.generate()), schedule="linear")
        a.run()
        b.run()
        print >>f, "File = " + filename + "; Schedule = geometric; Run #" + str(j+1) + "; Global Minimum PE = " + str(a.system_potential_energy())
        print >>f, "File = " + filename + "; Schedule = linear; Run #" + str(j+1) + "; Global Minimum PE = " + str(b.system_potential_energy())
f.close()
