from common import Particle, System

# default parameters
TARGET_TEMP = 1
RHO = 0.01
TOLERANCE = 0.001

# define optimization
class LocalMin_LJ(System):
    def __init__(self, filename, toler=TOLERANCE, temp=TARGET_TEMP, density=RHO):
        System.__init__(self, filename, temp, density)
        
        # set local optimization parameters
        self.LJ_lambda = 0.1
        self.LJ_energy_tolerance = toler
        self.energy_change = 1.0

    def run(self):
        while self.energy_change >= self.LJ_energy_tolerance:
            # save old positions and energies
            old_positions_table = self.collect_positions()
            old_energies_table = self.collect_energies()
            
            # move particles
            for particle in self.particles:
                (fx, fy, fz) = particle.force()
                particle.add_position(self.LJ_lambda*fx, self.LJ_lambda*fy, self.LJ_lambda*fz)

            # update and calculate new energies
            self._update_energies()
            old_energy, new_energy = sum(old_energies_table)/2.0, sum(self.collect_energies())/2.0
            # print str(old_energy) +"  "+str(new_energy)
            
            # decide to keep new configuration
            if new_energy < old_energy:
                self.LJ_lambda *= 1.2
                # for efficiency, forces will be updated ONLY if move is accepted
                self._update_forces()
                self.energy_change = abs(new_energy - old_energy)
            else:
                self.LJ_lambda *= 0.5
                self._restore_old_positions(old_positions_table)
                self._restore_old_energies(old_energies_table)

###########################
# RUN CODE
system_energies = []
f = open('hw4.1.results','w')

for i in range(5):
    filename = 'coord.' + str(i+1)
    a = LocalMin_LJ(filename)
    a.run()
    system_energies.append(a.system_potential_energy())
    print >>f, "File = " + filename + "; Local Minimum PE = " + str(system_energies[i])
print >>f, "Global Minimum PE = " + str(min(system_energies))
f.close()
