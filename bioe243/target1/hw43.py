from math import exp, sqrt
from common import Ran3, Lattice

# default parameters
SEED = 845340327839343
INIT_TEMP = 1.5

# define optimization
class SimulatedAnneal_Lattice(Lattice):
    def __init__(self, filename='/home/stuart/Desktop/BE243/Project/Target 1/Lprotein.36', temp=INIT_TEMP, r_seed=SEED):
        Lattice.__init__(self, filename, temp)
        self.rand2 = Ran3(r_seed)
        
    # one MC step (self.N MC moves)
    def MC_step():
        nbr_of_steps = 0
        for index in range(self.N):
            # get the mc move
            (coords, move_type) = self.mc_move(index)
            if coords != []:
                # if a move is available, save the old system
                old_sys_table = self.save_old_system()

                # make the move
                self._move_bead_to(index, coords[0])
                if move_type == "crankshaft":
                    self._move_bead_to(index+1, coords[1])

                self._update_energies()

                if not self.mc_accept(sum(old_sys_table[1]), sum(self.collect_energies())):
                    self._restore_old_system(old_sys_table)
                if self.mc_accept(sum(old_sys_table[1]), sum(self.collect_energies())):
                    nbr_of_steps += 1    
                    
                    
    
    def run(self):
        while self.temperature > 0.001:
            old_system_table = self.save_old_system()

            # Run MC steps
            steps = int(10 / self.temperature) #???
            for k in range(steps):
                self.MC_step()
            
            # averages - choose new config over old
            energy_change = self.system_potential_energy() - sum(old_system_table[1])/2.0
            if energy_change > 0.0 and self.rand2.generate() >= exp(-energy_change/self.temperature):
                # Revert back to state prior to simulating at current temperature,
                # since the new state has a higher energy AND rand is NOT < exp(-dE)
                self._restore_old_system(old_system_table)
            
            # run cooling schedule
            self.temperature *= 0.99
            
            print self.system_potential_energy()

###########################
# RUN CODE

#~ f = open('hw4.3.results','w')
seed = Ran3(1249834071)
self = []
V = SimulatedAnneal_Lattice.MC_step()
print V
#~ for i in range(3):
#~ a = SimulatedAnneal_Lattice(r_seed=int(1E12*seed.generate()))
#~ a.run()
    #~ print >>f, "File = " + 'protein_lattice_config_1.txt' + "; Schedule = geometric; Run #" + str(i+1) + "; Global Minimum PE = " + str(a.system_potential_energy())
#~ f.close()
