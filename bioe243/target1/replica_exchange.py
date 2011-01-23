#<<<<<<< Updated upstream
from random import random, randint, randrange
from math import exp
from copy import deepcopy # DEEPCOPY copies the entire memory space of an object and its sub-attributes
from energy_function import energy_function, convert_sequence
from heat_capacity import heat_capacity_function, MC_step

# assumes the temperatures are arranged in order from lowest to highest
def replica_exchange(configurations, temperatures, sequence, number_of_cycles):

    number_of_temperature_configurations = len(temperatures)
    
    # Initialize array_of_V's and record initial energies at each temperature
    # ex. [[V1, V2], [V3, V4]]
    array_of_Vs = []
    for i in xrange(number_of_temperature_configurations):
        array_of_Vs.append([energy_function(configurations[i], sequence)])
    
    for i in xrange(number_of_cycles):
        print "Running Cycle " + str(i)
        
        # Propose a replica exchange half of the time
        if random() < 0.5:
            '''
            # select ANY two different temperatures
            k = j = randrange(0, number_of_temperature_configurations)
            while j == k:
                j = randrange(0, number_of_temperature_configurations)
            '''    
            # Select two different temperatures ADJACENT to each other
            k = j = randrange(0, number_of_temperature_configurations)
            if k == 0:
                j += 1
            elif k == number_of_temperature_configurations-1 or random() < 0.5:
                j -= 1
            else:
                j += 1
            
            delta_beta = 1/temperatures[k]-1/temperatures[j];
            delta_energy = array_of_Vs[k][i-1] - array_of_Vs[j][i-1];
        
            # If we accept with acceptance min(1, exp(dU*dB)), then we exchange the systems (assumes that in C, it just involves pointer dereference swapping)
            if random() < min(1, exp(delta_beta*delta_energy)):
                #temperatures[k], temperatures[j] = temperatures[j], temperatures[k]
                configurations[k], configurations[j] = configurations[j], configurations[k]

        print energy_function(configurations[0], sequence)

        # For each temperature configuration
        for w in range(number_of_temperature_configurations):
            # Perform MC step on that configuration
            configurations[w],energy = MC_step(configurations[w], sequence, temperatures[w])

            # Save energy for that configuration at that time step
            array_of_Vs[w].append(energy)

    # when done, return the array_of_Vs
    return array_of_Vs
            
if __name__ == '__main__':
    filename = 'protein.36'
    #{
    import sys
    if len(sys.argv)>1:
        filename = sys.argv[1]
    #}
    configuration = []
    r_file = open(filename, 'r')
    r_line = r_file.readline()
    while r_line != '':
        r = tuple(int(p) for p in r_line.strip().replace(',', ' ').split())
        #r = Point(coord[0],coord[1],coord[2])
        configuration.append(r)
        #~ print configuration
        r_line = r_file.readline()
    r_file.close()

    # set up temperatures
    temperature = 0.10
    temperatures = []
    while temperature <= 1.1:
        temperatures.append(temperature)
        temperature += 0.05

    # set up configuration copies, all same initially
    configurations = []
    while len(configurations) < len(temperatures):
        configurations.append(list(configuration))

    # set up sequence and number_of_cycles
    sequence = convert_sequence('HHPHHPPHPPPHPPHPHHHHPPPPHHPPHPHHHHPP')
    number_of_cycles = 100 # CAN INCREASE THIS LATER
    
    # call replica exchange
    array_of_Vs = replica_exchange(configurations, temperatures, sequence, number_of_cycles)

    # print out the final energies (the last energy on array) and heat capacity to replica_exchange.py
    f = open('replica_exchange.results','w')
    print >>f, "HEAT CAPACITIES - REPLICA EXCHANGE IS USED \nNumber of Cycles = " + str(number_of_cycles)
    for c in range(len(configurations)):
        #print "Configuration #" + str(c+1) + \
        #      "; Temp = " + temperatures[c] + \
        #      "; Energy = " + str(array_of_Vs[c][number_of_cycles]) + \
        #      "; Heat Capacity = " + str(heat_capacity_function(array_of_Vs[c][number_of_cycles], number_of_cycles, temperatures[c]))
        print >>f, "Configuration %d; Temp = %f; Energy = %f; Heat Capacity = %f" % (
                    c+1,
                    temperatures[c],
                    array_of_Vs[c][-1],
                    heat_capacity_function(array_of_Vs[c], number_of_cycles, temperatures[c])
                )
    f.close()
