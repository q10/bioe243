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
                    direction_code = str(int(self.rand2.generate()*6))
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

    def evaluate_fitness(self):
        self._update_energies()

# define optimization algorithm
class GeneticAlgorithm:
    def __init__(self, filename='protein_lattice_config_1.txt', temp=INIT_TEMP, r_seed=SEED):
        self.population = []
        
        self.rand3 = Ran3(r_seed)
        self.crossover_probability = 0.5
        self.mutation_probability = 0.1

        self.construct_population()

        # calculate fitness levels of all members of population
        for structure in self.population:
            structure.evaluate_fitness()

        # for debugging
        for i in range(len(self.population)):
            print self.population[i].system_potential_energy()

    def construct_population(self):
        for i in range(200):
            self.population.append(Structure(seed=self.rand3.generate()))
        self.run_natural_selection()
        self.population = self.population[:100] # select 100 to begin with
    
    def run_natural_selection(self):
        self.population = filter(self.filter_out_different_length_structure, self.population)
        
    def filter_out_different_length_structure(self, structure):
        return len(structure.beads) == 27

    def individual_fitness(self, index):
        return self.population[index].system_potential_energy()

    def population_total_fitness(self):
        individual_fitnesses = []
        for i in range(len(self.population)):
            individual_fitnesses.append(self.individual_fitness(i))
        return sum(individual_fitnesses)
    
    def population_average_fitness(self):
        return self.population_total_fitness()/len(self.population)
                    
    def single_crossover(self, index_j, index_k):
        v1, v2 = self.population[index_j].vector, self.population[index_k].vector
        (v1, v2) = self.single_crossover_string(v1, v2)
        self.population[index_j].vector, self.population[index_k].vector = v1, v2

    def double_crossover(self, index_j, index_k):
        v1, v2 = self.population[index_j].vector, self.population[index_k].vector
        (v1, v2) = self.double_crossover_string(v1, v2)
        self.population[index_j].vector, self.population[index_k].vector = v1, v2

    def single_crossover_string(self, v1, v2):
        i = int(self.rand3.generate()*27)
        A1, B1, A2, B2 = v1[:i], v2[:i], v1[i:], v2[i:]
        return (A1+B2, B1+A2)

    def double_crossover_string(self, v1, v2):
        i, j = int(self.rand3.generate()*27), int(self.rand3.generate()*27)
        if j < i:
            i,j = j,i
        A1, B1, A2, B2, A3, B3 = v1[:i], v2[:i], v1[i:j], v2[i:j], v1[j:], v2[j:]
        return (A1+B2+A3, B1+A2+B3)

    def select_two_parents(self, best, biggest_energy_difference, total_fitness):
        # based probability on the fitness function in page 27 of Lecture 14 slides
        '''
        i = int(self.rand3.generate()*len(self.population))
        inverse_prob = abs(self.individual_fitness(i)-self.individual_fitness(best))/biggest_energy_difference
        while self.rand3.generate() >= 1 - inverse_prob:
            i = int(self.rand3.generate()*len(self.population))
            inverse_prob = abs(self.individual_fitness(i)-self.individual_fitness(best))/biggest_energy_difference
            
        j = int(self.rand3.generate()*len(self.population))
        inverse_prob = abs(self.individual_fitness(j)-self.individual_fitness(best))/biggest_energy_difference
        while j != i and self.rand3.generate() >= 1 - inverse_prob:
            j = int(self.rand3.generate()*len(self.population))
            inverse_prob = abs(self.individual_fitness(j)-self.individual_fitness(best))/biggest_energy_difference
        return (i, j)
        '''
        return (int(self.rand3.generate()*len(self.population)), int(self.rand3.generate()*len(self.population)))

    def apply_crossover_and_create_children(self, index_j, index_k):
        if self.rand3.generate() < self.crossover_probability:
            if self.rand3.generate() < 0.7:
                self.single_crossover(index_j, index_k)
            else:
                self.double_crossover(index_j, index_k)
        # the above methds crosses over the vectors only; no need to create new objects
        # in either case, we return both structures back, crossed or not
        return [self.population[index_j], self.population[index_k]]
    
    def apply_mutation(self, children):
        for child in children:
            if self.rand3.generate() < self.mutation_probability:
                # run point mutation
                point = int(self.rand3.generate()*27)
                A, B, C = child.vector[:point], str(int(self.rand3.generate()*6)), child.vector[:point+1]
                child.vector = A + B + C
        return children
                
    def find_worst_and_best_indices(self):
        worst, best = 0, 0
        for i in range(1, len(self.population)):
            tmp = self.individual_fitness(i)
            if tmp < self.individual_fitness(best):
                best = i
            if tmp > self.individual_fitness(worst):
                worst = i
        return (worst, best)
    
    def run(self):
        old_fitness_score, new_fitness_score = 0, 100
        while abs(old_fitness_score - new_fitness_score) > 0.01:

            old_fitness_score = new_fitness_score
            
            # initialize new population        
            new_population = []
            
            (worst, best) = self.find_worst_and_best_indices()
            biggest_energy_difference = abs(self.individual_fitness(worst) - self.individual_fitness(best))
            total_fitness = self.population_total_fitness()

            # save the best structure into the next generation (Elitism)
            new_population.append(self.population[best])

            while len(new_population) < len(self.population):

                # select parents based on fitness scores
                (j, k) = self.select_two_parents(best, biggest_energy_difference, total_fitness)

                # with probability, create children
                children = self.apply_crossover_and_create_children(j, k)

                # with probability, mutate them
                children = self.apply_mutation(children)

                # regrow chain, helps determine which children are viable
                for child in children:
                    child.regenerate_lattice_from_vector()

                # filter chilren for viability (ex. structures have overlapping beads, etc)
                children = filter(self.filter_out_different_length_structure, children)

                # place children in new population
                new_population.extend(children)

            # ensures we get new population of same size
            self.population = new_population[:100]

            # re-calculate fitness levels of all members of population
            for structure in self.population:
                structure.evaluate_fitness()

            new_fitness_score = self.population_average_fitness()
            #print self.population[0].system_potential_energy()
            # for debugging
            print 'running'
            print self.individual_fitness(best)
            '''for i in range(len(self.population)):
                print self.population[i].system_potential_energy()'''


###########################
# RUN CODE
a=GeneticAlgorithm()
a.run()
'''
f = open('hw4.4.results','w')
seed = Ran3(1249834071)

for i in range(3):
    a = GeneticAlgorithm(r_seed=int(1E12*seed.generate()))
    a.run()
    print >>f, "File = " + 'protein_lattice_config_1.txt' + "; Genetic Algorithm; Run #" + str(i+1) + "; Global Minimum PE = " + str(a.system_potential_energy())
f.close()
'''
