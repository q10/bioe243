from math import exp, sqrt
from common import Ran3, Lattice

# default parameters
SEED = 845340327839343
INIT_TEMP = 1.5

class Structure(Lattice):
    def __init__(self, v, filename='protein_lattice_config_1.txt'):
        self.vector = v # 27 chars long
        self.N = len(self.vector)
        
        self.beads = []
        self.__create_beads()

        self._update_energy()

    def read_in_bead_type(self, filename):
        bead_type = []
        r_file = open(filename, 'r')
        r_line = r_file.readline()
        while r_line != '':
            r = r_line.strip().replace(',', ' ').split()
            bead_type.append(r[0])
            r_line = r_file.readline()
        r_file.close()
        return bead_type
        
    def __create_beads(self, filename):
        bead_type = self.read_in_bead_type(filename)
        # start coords at origin
        x, y, z = 0, 0, 0
        for i in range(self.N):
            (x, y, z) = self.next_position(x, y, z, self.vector[i])
            self.beads.append(AABead(x, y, z, bead_type[i]))
        
    def regenerate_lattice_from_vector(self):
        x, y, z = 0, 0, 0
        for i in range(self.N):
            (x, y, z) = self.next_position(x, y, z, self.vector[i])
            self.beads[i].set_position(x, y, z)

    def next_position(self, x, y, z, direction):
        if direction == '0':
            y += 1
        if direction == '1':
            y -= 1
        if direction == '2':
            x -= 1
        if direction == '3':
            x += 1
        if direction == '4':
            z -= 1
        if direction == '5':
            z += 1
        return (x, y, z)

# define optimization
class GeneticAlgorithm:
    def __init__(self, filename='protein_lattice_config_1.txt', temp=INIT_TEMP, r_seed=SEED):
        self.population = []
        self.population.append(Structure(temp_vector))
        self.rand2 = Ran3(r_seed)
    def individual_fitness(self):
    def population_average_fitness(self):
    def mate_pairs(self):
    def apply_crossover(self, index_j, index_k):
        v1, v2 = self.population[index_j].vector, self.population[index_k].vector
        (v1, v2) = self.crossover_string(v1, v2, int(rand2.generate()*27)]
        self.population[index_j].vector, self.population[index_k].vector = v1, v2
        
    def apply_mutation(self, index):
        self.population[index].vector[int(rand2.generate()*27)] = str(int(rand2.generate()*5))

    def crossover_string(self, v1, v2, v_index):
        first, second, tmp, tmp2 = v1[:v_index], v2[:v_index], v1[v_index:], v2[v_index:]
        return (first+tmp2, second+tmp)
        
    def run(self):
        while self.temperature > 0.001:
            old_system_table = self.save_old_system()

            # Run MC steps
            steps = int(10 / self.temperature)
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
            
            # print self.system_potential_energy()

###########################
# RUN CODE

f = open('hw4.4.results','w')
seed = Ran3(1249834071)

for i in range(3):
    a = GeneticAlgorithm_Lattice(r_seed=int(1E12*seed.generate()))
    a.run()
    print >>f, "File = " + 'protein_lattice_config_1.txt' + "; Schedule = geometric; Run #" + str(i+1) + "; Global Minimum PE = " + str(a.system_potential_energy())
f.close()
