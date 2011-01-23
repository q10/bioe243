#################
## Monte Carlo ##
#################

from random import random
from math import exp
from energy_function import energy_function
#### Note: I have not modified this code for this target, but I can confirm it works as-is
## it expects my previous naming, such as calc_energy instead of energy_function, etc

# sorry, more weird lambdas. I've included the explanation in this file for nbr_gen_3d

nbr_gen_3d = lambda i,p,nbr: p + ((1-((nbr&1)<<1)) if (((nbr&6)>>1)==i) else 0)
get_neighbors = lambda pos: [tuple(nbr_gen_3d(i,p,nbr) for i,p in enumerate(pos)) for nbr in range(6)]
is_empty = lambda pos,config: pos not in config
get_empties = lambda nbrs,config: [nbr for nbr in nbrs if is_empty(nbr,config)]
is_bend = lambda before,after: len([ 1 for i in range(3) if abs(before[i]-after[i])==1 ])==2

def find_2(f, seq, seq2):
    """Return first index in sequence where f(item,item2) == True."""
    for i in xrange(len(seq)):
        if f(seq[i],seq2[i]):
            return i

def find_flip(cur,before,after):
    flat_direction = find_2(lambda a,b:a==b,after,before)
    bend_direction = -1 # find_2(lambda c,a:c==a, cur ,after)
    for i in range(3):
        if i==flat_direction: continue
        if cur[i]==after[i]: bend_direction = i
    new_pos = []
    for i in range(3):
        if i==flat_direction:
            new_pos.append(cur[i])
        elif i==bend_direction:
            new_pos.append(before[i])
        else:
            new_pos.append(after[i])
    return tuple(new_pos)

def get_crank_moves(residue, configuration, before, after, flip):
    # get other crankshaft member
    if flip == configuration[residue-2]: #before us
        #if sequential: None # we've already checked this one out
        other = residue-1
        rotators = after,flip
    elif flip == configuration[residue+2]:
        other = residue+1
        rotators = before,flip
    else:
        return []
    # check if crank-to spaces are open
    all_nbrs = [get_neighbors(rotators[0]),get_neighbors(rotators[1])]
    possible_moves = [(nbr1,nbr2) for nbr1,nbr2 in zip(all_nbrs[0],all_nbrs[1]) if is_empty(nbr1,configuration) and is_empty(nbr2,configuration)]

    if not len(possible_moves): return []
    else:
        crank_moves = [((residue,other),move) for move in possible_moves]
        return crank_moves

# returns a tuple(list(tuple(tuple(),tuple())),string)
# which means tuple(moves,meaning) where moves is a list of tuple(residues,targets)
# to accept a move, all residues must be moved to respective targets
def find_move(residue, configuration, sequence_len, sequential=True):
    # (cheapest) check for end move
    if residue==0 or residue==sequence_len-1:
        return get_end_moves(residue, configuration)
    # check for corner flip, we are guaranteed a before & after
    cur_pos = configuration[residue]
    before,after = configuration[residue-1],configuration[residue+1]
    moves = []
    if is_bend(before,after): # this is a bend
        flip = find_flip(cur_pos,before,after)
        if not flip in configuration:
            moves = [((residue,),(flip,))]
        # lastly, check for crankshaft
        # can only crankshaft if we're in a bend with no corner flips possible
        # and if we're not toward an end
        if residue==1 or residue==sequence_len-2: return moves

        moves.extend(get_crank_moves(residue, configuration, before, after, flip))
    return moves


def calc_acceptance(curr_E,next_E,sequence,T_star):
    if next_E<=curr_E:
        return True
    deltaV = next_E-curr_E
    factor = exp(-deltaV/T_star)
    if factor>random():
        return True
    return False

def get_end_moves(residue, configuration):
    if residue==0:  p_residue = residue+1
    else:           p_residue = residue-1
    cur_pos = configuration[p_residue]
    empties = get_empties(get_neighbors(cur_pos),configuration)
    end_moves = [((residue,),(e,)) for e in empties]
    return end_moves


#def mmc(ncycles=1):
if __name__=='__main__':
    sequence = (0,1,0,0,1,1,1,1,0,0,0,0,1)
    sequence_len = len(sequence)
    current_configuration = [(0,0,0),(1,1,1),(2,2,2)]
    residue = 1
    configuration = current_configuration

    cur_pos,before, after = configuration[residue],configuration[residue-1],configuration[residue+1]
    flip = find_flip(cur_pos,before, after)
    moves = get_crank_moves(residue, configuration, before, after, flip)

    moves = find_move(residue,current_configuration,sequence_len)


    if calc_acceptance(1.0, 2.0, sequence, 2.0):
        print moves

