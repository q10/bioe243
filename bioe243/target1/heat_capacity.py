from energy_function import energy_function, convert_sequence
from montecarlo import calc_acceptance, apply_move, find_move

def MC_step(present_configuration, sequence, T_star):
    # run over all beads in configuration
    present_V = energy_function(present_configuration,sequence)
    for m in range(len(sequence)):
        moves = find_move(m, present_configuration, len(sequence))
        if not moves: continue

        next_configuration = apply_move(moves, present_configuration)
        next_V = energy_function(next_configuration, sequence)
        if calc_acceptance(present_V, next_V, sequence, T_star):
            present_configuration = next_configuration
            present_V = next_V

    return present_configuration, present_V

# V = an array of energy of system after each MC step (i.e. V[i] will be configuration energy after the ith MC step)
# the MC code is responsible for saving values into V after each MC step
def heat_capacity_function(V, number_of_cycles, temperature):
    V2 = [] # <V^2>
    for energy in V:
        V2.append(energy**2)
    V, V2 = sum(V)/number_of_cycles, sum(V2)/number_of_cycles
    return (V2 - (V**2)) / (temperature**2)

# calculates heat capacity of the configurations at different temperatures - NOT USED, MAY BE OBSOLETED
def heat_capacities_of_all_temperature_configurations(array_of_Vs, number_of_cycles, temperatures):
    heat_capacities = []
    for i in range(temperatures):
        heat_capacities.append(heat_capacity_function(array_of_Vs[i], number_of_cycles, temperatures[i]))
    return heat_capacities

# function essentially the same as replica_exchange, but w/o the replica exchange
def run_multiple_MC(configurations, temperatures, sequence, number_of_cycles):
    number_of_temperature_configurations = len(temperatures)
    
    # Initialize array_of_V's and record initial energies at each temperature
    # ex. [[V1, V2], [V3, V4]]
    array_of_Vs = []
    for i in range(number_of_temperature_configurations):
        array_of_Vs.append([energy_function(configurations[i], sequence)])
    
    for i in range(number_of_cycles):
        print "running cycle " + str(i)
        # For each temperature configuration
        for w in range(number_of_temperature_configurations):

            # Perform MC step on that configuration
            configurations[w],energy = MC_step(configurations[w], sequence, temperatures[w])

            # Save energy for that configuration at that time step
            array_of_Vs[w].append(energy)


    # when done, return the array_of_Vs
    #~ print "Got V of %s" % array_of_Vs
    return array_of_Vs

# remove this and use MC_step() from replica_exchange.py
'''
def MC_step(configuration, sequence, temperature):
    possible_moves = []
    for r in xrange(len(sequence)):
        moves = find_move(r,configuration,len(sequence))
        if not moves: continue
        # if we need to try a move at each residue, use these 5 lines
        #next_config = apply_move(moves,configurations[w])
        ## accept?
        #if calc_acceptance(configurations[w],next_config,temperatures[w]):
        #    # accept!
        #    configurations[w] = next_config
        
        # if we try only one move along the whole residue, use these lines:
        possible_moves.extend(moves)
    next_config = apply_move(possible_moves,configuration)
    if calc_acceptance(configuration,next_config,sequence, temperature):
        configuration = next_config
        #~ print "Accepted, passed through"
        #~ print configuration
    return configuration
'''        
if __name__=='__main__':
    filename = 'protein.36'
    #{
    #~ import sys
    #~ if len(sys.argv)>1:
        #~ filename = sys.argv[1]
    #}
    configuration = []
    r_file = open(filename, 'r')
    r_line = r_file.readline()
    while r_line != '':
        r = tuple(map(int, r_line.strip().replace(',', ' ').split()))
        #r = Point(coord[0],coord[1],coord[2])
        configuration.append(r)
        #~ print configuration
        r_line = r_file.readline()
    r_file.close()

    # set up temperatures
    temperature = 0.10
    temperatures = []
    while temperature <= 1.105:
        temperatures.append(temperature)
        temperature += 0.05
    print temperatures

    # set up configuration copies, all same initially
    configurations = []
    while len(configurations) < len(temperatures):
        configurations.append(list(configuration))

    # set up sequence and number_of_cycles
    #sequence = convert_sequence('HHPHHPPHPPPHPPHPHHHHPPPPHHPPHPHHHHPP')
    sequence = convert_sequence('HHHPPHHPHHHHPHHPPHPHHPHPPHPPPHHPPPPP')
    number_of_cycles = 10000 # CAN INCREASE THIS LATER

    # call MC cycle
    array_of_Vs = run_multiple_MC(configurations, temperatures, sequence, number_of_cycles)

    # print out the final energies (the last energy on array) and heat capacity to heat_capacity.results
    f = open('heat_capacity.results','w')
    print >>f, "HEAT CAPACITIES - NO REPLICA EXCHANGE HERE \nNumber of Cycles = " + str(number_of_cycles)
    for c in range(len(configurations)):
        #print "Configuration #" + str(c+1) + \
        #      "; Temp = " + str(temperatures[c]) + \
        #      "; Energy = " + str(array_of_Vs[c][number_of_cycles]) + \
        #      "; Heat Capacity = " + str(heat_capacity_function(array_of_Vs[c][number_of_cycles], number_of_cycles, temperatures[c]))
        print >>f, "Configuration #%d; Temp = %.2f; Energy = %f; Heat Capacity = %f" % (
                    c+1,temperatures[c],
                    array_of_Vs[c][-1],
                    heat_capacity_function(array_of_Vs[c], number_of_cycles, temperatures[c])
            )
    f.close()
