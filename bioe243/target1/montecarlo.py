#################
## Monte Carlo ##
#################

from random import choice
from energy_function import energy_function
from montecarlo_ss import find_move,calc_acceptance

def apply_move(moves, config):
    move = choice(moves)
    new_config = list(config)
    for i in range(len(move[0])):
        new_config[move[0][i]]=move[1][i]
    return new_config

#def mmc(ncycles=1):
if __name__=='__main__':
    sequence = (0,1,0,0,1,1,1,1,0,0,0,0,1)
    sequence_len = len(sequence)
    current_configuration = [(0,0,0),(1,1,1)]
    residue = 1

    moves = find_move(residue,current_configuration,sequence_len)

    current_configuration = apply_move(moves,current_configuration)
    if calc_acceptance(1.0, 2.0, sequence, 2.0):
        print moves

