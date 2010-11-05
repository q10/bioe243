from math import exp, sqrt
from common import Ran3, AABead, Lattice

# default parameters
SEED = 845340327839343
INIT_TEMP = 1.5

class Structure(Lattice):
    def __init__(self, filename='protein_lattice_config_1.txt', seed=SEED):
        self.vector = "" # 27 chars long
        self.beads = []

        self.rand2 = Ran3(seed)
        self.__create_beads(filename)

        self._update_energies()
        
    def __create_beads(self, filename):
        r_file = open(filename, 'r')
        r_line = r_file.readline()

        # start coords at origin
        current_coords = (0, 0, 0)

        # grow chain, and terminate early if chain end is "stuck" in a spot with no free adjacent positions
        while r_line != '':            
            if self.spaces_left_for_next_bead(current_coords) == None:
                break
            else:
                next_coords = current_coords
                direction_code = ''
                while self.find_bead_at_this_coordinate(next_coords) != None:
                    direction_code = str(int(self.rand2.generate()*5))
                    next_coords = self.next_position(current_coords, direction_code)
            
                bead_typ = r_line.strip().replace(',', ' ').split()[0]
                self.beads.append(AABead(next_coords[0], next_coords[1], next_coords[2], bead_typ))
                current_coords = next_coords
                self.vector += direction_code
                r_line = r_file.readline()
        r_file.close()

    def spaces_left_for_next_bead(self, current_coords):
        neighbor_coords = self.neighbor_coords_of_bead(AABead(current_coords[0], current_coords[1], current_coords[2], 'H'))
        for neighbor in neighbor_coords:
            if self.find_bead_at_this_coordinate(neighbor) == None:
                return True
        return None
        
    def regenerate_lattice_from_vector(self):
        current_coords = (0, 0, 0)
        for i in range(len(self.vector)):
            current_coords = self.next_position(current_coords, self.vector[i])
            if self.find_bead_at_this_coordinate(current_coords) == None:
                self.beads[i].set_position(current_coords[0], current_coords[1], current_coords[2])
            else:
                # throw errors by shortening string, genetic algorithm will dispose of this as "unfit," since it's less than 27 beads long
                self.vector = self.vector[:i]
                break

    def next_position(self, coords, direction):
        (x, y, z) = coords
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

# define optimization algorithm
class GeneticAlgorithm:
    def __init__(self, filename='protein_lattice_config_1.txt', temp=INIT_TEMP, r_seed=SEED):
        self.population = []
        self.rand3 = Ran3(r_seed)
        
        for i in range(100):
            self.population.append(Structure(seed=self.rand3.generate()))
        

    def individual_fitness(self, index):
        return self.population[index].system_potential_energy()
    
    def population_average_fitness(self):
        individual_fitnesses = []
        for i in range(len(self.population)):
            individual_fitnesses.append(self.individual_fitness(i))
        return sum(individual_fitnesses)/len(self.population)
        
    def mate_pairs(self):
    def apply_crossover(self, index_j, index_k):

    def single_crossover(self, index_j, index_k):
        v1, v2 = self.population[index_j].vector, self.population[index_k].vector
        (v1, v2) = self.single_crossover_string(v1, v2)
        self.population[index_j].vector, self.population[index_k].vector = v1, v2

    def double_crossover(self, index_j, index_k):
        v1, v2 = self.population[index_j].vector, self.population[index_k].vector
        (v1, v2) = self.double_crossover_string(v1, v2, int(rand2.generate()*27))
        self.population[index_j].vector, self.population[index_k].vector = v1, v2

    def apply_mutation(self, index):
        self.population[index].vector[int(rand2.generate()*27)] = str(int(rand2.generate()*5))

    def single_crossover_string(self, v1, v2):
        i = int(rand2.generate()*27)
        A1, B1, A2, B2 = v1[:i], v2[:i], v1[i:], v2[i:]
        return (A1+B2, B1+A2)

    def double_crossover_string(self, v1, v2):
        i, j = int(rand2.generate()*27), int(rand2.generate()*27)
        if j < i:
            i,j = j,i
        A1, B1, A2, B2, A3, B3 = v1[:i], v2[:i], v1[i:j], v2[i:j], v1[j:], v2[j:]
        return (A1+B2+A3, B1+A2+B3)
        
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
'''
###########################
# RUN CODE

f = open('hw4.4.results','w')
seed = Ran3(1249834071)

for i in range(3):
    a = GeneticAlgorithm(r_seed=int(1E12*seed.generate()))
    a.run()
    print >>f, "File = " + 'protein_lattice_config_1.txt' + "; Schedule = geometric; Run #" + str(i+1) + "; Global Minimum PE = " + str(a.system_potential_energy())
f.close()
'''
