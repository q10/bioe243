from energy_function import energy_function, convert_sequence
from montecarlo import calc_acceptance, apply_move, find_move
from random import random, randint, randrange
from math import exp
from copy import deepcopy # DEEPCOPY copies the entire memory space of an object and its sub-attributes

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
    return (present_configuration, sequence_traversals, nbr_moves)

if __name__=='__main__':
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
    start_configuration = configuration
    
    # set up temperature, target folding energy
    temperature = 0.65
    Vmin = -34.5

    # set up sequence and number_of_cycles
    #sequence = convert_sequence('HHPHHPPHPPPHPPHPHHHHPPPPHHPPHPHHHHPP')
    sequence = convert_sequence('HHHPPHHPHHHHPHHPPHPHHPHPPHPPPHHPPPPP')
    Max_number_of_cycles = 100 # CAN INCREASE THIS LATER, not currently used
    
    # counters for MPFT
    sequence_traversals = nbr_moves = 0
    trajectories = 3
    running_total = 0
    
    MFPT = []
    all_sequence_traversals = []
    
    V = energy_function(configuration, sequence)
    print 'Starting V = ' + str(V)
    for i in xrange(trajectories):
        count = 0
        while V > Vmin:
            count += 1
            #~ configurations[w] = MC_step(configurations[w], sequence, temperatures[w])
            #~ Result = (configuration, times down sequence, number of moves)
            Result = MC_step(configuration, sequence, temperature, sequence_traversals, nbr_moves)
            configuration = Result[0]
            V = energy_function(configuration, sequence)
            running_total = running_total + Result[2]
            if count % 1000 == 0:
                print '%f steps taken so far, at energy %f' % (running_total, V)
            #~ print running_total
        #~ print nbr_of_steps
        print "reached energy after %f steps" % (running_total)
        MFPT.append(running_total)
        #~ all_sequence_traversals.append(Result[1])
        sequence_traversals = nbr_moves = running_total = 0 # Reset counters for next trajectory
        configuration = start_configuration
        V = energy_function(configuration, sequence)
    print MFPT
    #~ print 'sequceT' +str(all_sequence_traversals)
    averageMFPT = ( sum(MFPT) / trajectories )
    print averageMFPT
    #~ for c in range(len(configurations)):
    print "Configuration #" + str(1) + \
          "; Temp = " + str(temperature) + \
          "; Energy = " + str(V) + \
          "; Average number of MC Steps = " + str(averageMFPT)
