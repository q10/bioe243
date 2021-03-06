from math import exp, sqrt
from common import Ran3, AABead, Lattice
from copy import deepcopy

# default parameters
SEED = 845340327839343
INIT_TEMP = 1.5

class Structure(Lattice):
    def __init__(self, filename='protein_lattice_config_1.txt', seed=SEED):
        self.vector = "" # 27 chars long
        self.beads = []

        self.rand2 = Ran3(seed)
        self.__create_beads(filename)

    def __str__(self):
        return self.vector
    
    def __create_beads(self, filename):
        r_file = open(filename, 'r')
        r_line = r_file.readline()

        # start coords at ONE unit from origin
        current_coords = (0, 0, 0)
        direction_code = str(int(self.rand2.generate()*6))
        next_coords = self.next_position(current_coords, direction_code)
        bead_typ = r_line.strip().replace(',', ' ').split()
        self.beads.append(AABead(next_coords[0], next_coords[1], next_coords[2], bead_typ[0]))
        current_coords = next_coords
        self.vector += direction_code
        r_line = r_file.readline()

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
            
                bead_typ = r_line.strip().replace(',', ' ').split()
                self.beads.append(AABead(next_coords[0], next_coords[1], next_coords[2], bead_typ[0]))
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
        for bead in self.beads:
            bead.set_position(0,0,0)
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
    def __init__(self, filename='protein_lattice_config_1.txt', temp=INIT_TEMP, r_seed=SEED, population=100, run_times=10000):
        self.population = []
        
        self.rand3 = Ran3(r_seed)
        self.crossover_probability = 0.5
        self.mutation_probability = 0.1
        self.population_size = population
        self.run_times=run_times

        self.construct_population()

        # calculate fitness levels of all members of population
        for structure in self.population:
            structure.evaluate_fitness()

    def construct_population(self):
        while len(self.population) < self.population_size:
            self.population.append(Structure(seed=self.rand3.generate()))
            self.run_natural_selection()
    
    def run_natural_selection(self):
        self.population = filter(self.filter_out_different_length_structure, self.population)
        
    def filter_out_different_length_structure(self, structure):
        return len(structure.vector) == 27

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
        child1, child2 = deepcopy(self.population[index_j]), deepcopy(self.population[index_k])
        child1.vector, child2.vector = v1, v2
        return [child1, child2]

    def double_crossover(self, index_j, index_k):
        v1, v2 = self.population[index_j].vector, self.population[index_k].vector
        (v1, v2) = self.double_crossover_string(v1, v2)
        child1, child2 = deepcopy(self.population[index_j]), deepcopy(self.population[index_k])
        child1.vector, child2.vector = v1, v2
        return [child1, child2]

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
        i = int(self.rand3.generate()*len(self.population))
        prob = abs(self.individual_fitness(i) / total_fitness)
        while self.rand3.generate() >= prob:
            i = int(self.rand3.generate()*len(self.population))
            prob = abs(self.individual_fitness(i) / total_fitness)

        j = int(self.rand3.generate()*len(self.population))
        prob = abs(self.individual_fitness(j) / total_fitness)
        while self.rand3.generate() >= prob or j == i:
            j = int(self.rand3.generate()*len(self.population))
            prob = abs(self.individual_fitness(i) / total_fitness)
        return (i, j)

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
        
        return (int(self.rand3.generate()*len(self.population)), int(self.rand3.generate()*len(self.population)))
        '''
    def apply_crossover_and_create_children(self, index_j, index_k):
        children = []
        if self.rand3.generate() < self.crossover_probability:
            if self.rand3.generate() < 0.7:
                children = self.single_crossover(index_j, index_k)
            else:
                children = self.double_crossover(index_j, index_k)
        else:
            # else return parents as members of the new population
            # deepcopy is an imported function that creates a copy of an object in memory
            children = [deepcopy(self.population[index_j]), deepcopy(self.population[index_k])]
        return children
    
    def apply_mutation(self, children):
        for child in children:
            if self.rand3.generate() < self.mutation_probability:
                # run point mutation
                point = int(self.rand3.generate()*27)
                A, B, C = child.vector[:point], str(int(self.rand3.generate()*6)), child.vector[point+1:]
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

    def find_second_best_index(self):
        best = 0
        for i in range(1, len(self.population)):
            tmp = self.individual_fitness(i)
            if tmp < self.individual_fitness(best):
                best = i
        second_best = 0
        for j in range(1, len(self.population)):
            tmp = self.individual_fitness(j)
            if tmp < self.individual_fitness(second_best) and j != best:
                second_best = j
        return second_best

    def print_population(self):
        print "CURRENT POPULATION STRINGS"
        for i in range(len(self.population)):
            print str(self.population[i])
    
    def run(self):
        for t in range(self.run_times):
            # initialize new population and set up other variables for later use in loop
            new_population = []
            (worst, best) = self.find_worst_and_best_indices()
            biggest_energy_difference = abs(self.individual_fitness(worst) - self.individual_fitness(best))
            total_fitness = self.population_total_fitness()

            # save the two best structures into the next generation (Elitism)
            new_population.append(self.population[best])
            new_population.append(self.population[self.find_second_best_index()])
            
            while len(new_population) < len(self.population):
                '''print "new pop length: "+str(len(new_population))'''

                # select parents based on their fitness scores
                (j, k) = self.select_two_parents(best, biggest_energy_difference, total_fitness)

                #print "begin:\t\t" + self.population[j].vector + " " + self.population[k].vector

                # with crossover_probability, create children
                children = self.apply_crossover_and_create_children(j, k)
                # print "children:\t" + children[0].vector + " " + children[1].vector

                # with mutation_probability, mutate them
                children = self.apply_mutation(children)
                # print "children:\t" + children[0].vector + " " + children[1].vector

                # regrow chain, helps determine which children are viable
                for child in children:
                    child.regenerate_lattice_from_vector()

                '''
                for child in children:
                    print "children okay1 " + str(child)
                print "\n"
                '''

                # filter chilren for "viability" (ex. structures have overlapping beads, structure is stuck and cannot grow fully, etc)
                children = filter(self.filter_out_different_length_structure, children)

                '''
                for child in children:
                    print "children okay2 " + str(child)
                print "\n"
                '''

                # place children in new population
                new_population.extend(children)

            '''
            print "NEW_POP"
            for i in range(len(new_population)):
                print str(new_population[i])
            '''
            
            # ensures we get new population of same size
            self.population = new_population[:self.population_size]

            # re-calculate fitness levels of all members of population
            for structure in self.population:
                structure.evaluate_fitness()

            '''
            print self.population[0].system_potential_energy()
            (worst, best) = self.find_worst_and_best_indices()
            print 'running' + str(self.individual_fitness(best))
            for i in range(len(self.population)):
                print self.population[i].system_potential_energy()
            '''

###########################
# RUN CODE

f = open('hw4.4.results','w')
seed = Ran3(1249834071)

for w in range(3):
    a = GeneticAlgorithm(r_seed=int(1E12*seed.generate()), population=200*(w+1), run_times=500)
    a.run()
    (worst, best) = a.find_worst_and_best_indices()
    print >>f, "File = " + 'protein_lattice_config_1.txt' + "; Genetic Algorithm\nRun #" + str(w+1) + \
          "\nPopulation Size = " + str(len(a.population)) + \
          "\nGenerations Iterated = " + str(a.run_times) + \
          "\nAverage Fitness (PE) of Solution Population = " + str(a.population_average_fitness()) + \
          "\nBest Global Minimum PE = " + str(a.individual_fitness(best)) + "\n\n"
f.close()
