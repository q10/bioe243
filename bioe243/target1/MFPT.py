from energy_function import energy_function, convert_sequence
from genetic_lattice import encode_lattice
from montecarlo import calc_acceptance, apply_move, find_move
from random import random, randint, randrange
from math import exp
from copy import deepcopy # DEEPCOPY copies the entire memory space of an object and its sub-attributes
from global_minimum import print_at

def MC_step(present_configuration, sequence, T_star,sequence_traversals, nbr_moves):
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
            nbr_moves += 1
            #~ print 'Moved! Number of moves = %f and presentV = %f .' % (nbr_moves, present_V)
    sequence_traversals += 1
    return present_configuration, sequence_traversals, nbr_moves, present_V

if __name__=='__main__':
    filename = 'protein.36'
    import sys
    if len(sys.argv)>1:
        filename = sys.argv[1]

    if 'second' not in sys.argv:
        sequence = convert_sequence('HHPHHPPHPPPHPPHPHHHHPPPPHHPPHPHHHHPP')
    else:
        sequence = convert_sequence('HHHPPHHPHHHHPHHPPHPHHPHPPHPPPHHPPPPP')
    sequence_len = len(sequence)

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
    start_configuration = configuration
    
    # set up temperature, target folding energy
    temperature = 0.65
    Vmin = float(sys.argv[2]) if len(sys.argv)>2 else -29

    # set up sequence and number_of_cycles
    Max_number_of_cycles = 100 # CAN INCREASE THIS LATER, not currently used
    
    # counters for MPFT
    sequence_traversals = nbr_moves = 0
    trajectories = 3
    running_total = 0
    
    MFPT = []
    all_sequence_traversals = []
    
    V = energy_function(configuration, sequence)
    for i in xrange(trajectories):
        sequence_traversals = nbr_moves = 0 # Reset counters for next trajectory
        configuration = start_configuration
        V = energy_function(configuration, sequence)
        print_at('Starting V = %f' % V,1)
        print_at('Trajectory %d' % i,2)
        j = 0
        lowest_V = 0
        while V-0.0001 > Vmin: # floating point errors...
            j += 1
            print_at('Move %d' % j,3)
            print_at('With V = %f' % V,4)
            #~ configurations[w] = MC_step(configurations[w], sequence, temperatures[w])
            #~ Result = (configuration, times down sequence, number of moves)
            configuration, sequence_traversals, nbr_moves, V = MC_step(
                     configuration, sequence, temperature, sequence_traversals, nbr_moves)
            if V < lowest_V:
                lowest_V = V
                lowest_conf = configuration
                print_at('best config w/ E=%f: %s' % (lowest_V,str(encode_lattice(lowest_conf)).replace(' ','')),5)
            #~ print Result
        #~ print nbr_of_steps
        print_at("reached energy after %f steps with config %s" % (
                    nbr_moves,str(encode_lattice(configuration)).replace(' ','')
            ),7)
        MFPT.append(nbr_moves)
        all_sequence_traversals.append(sequence_traversals)
    print MFPT
    #~ print 'sequceT' +str(all_sequence_traversals)
    averageMFPT = ( sum(MFPT) / trajectories )
    print averageMFPT
    print "Configuration #%d; Temp = %f; Energy = %f; Average number of MC Steps = %d" % (
            1,
            temperature,
            V,
            averageMFPT
        )
