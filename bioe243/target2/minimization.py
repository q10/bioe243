def check_local_minimum(seq, moves, orig_config):
    #print "(local) checking... %s" % seq
    # 1k MC steps
    from minimization import apply_move, find_move
    from energy_function import energy_function
    from genetic_lattice import print_at
    init_energy = energy_function(orig_config, seq)
    config = orig_config
    for move in moves:
        nconfig = apply_move(move, config)
        if energy_function(nconfig, seq)<init_energy:
            #print_at("found lower energy with MC... %s" % (energy_function(nconfig, seq)),4)
            return (False, seq)
    #print_at("quick MC was a'ight %s" % seq,1)
    from minimization import calc_acceptance, apply_some_move
    sequence_len = len(seq)
    sequence_range = xrange(sequence_len)
    cE = energy_function(config, seq)
    #print_at("MC step?",2)
    for i in xrange(10**3):
        for res in sequence_range:
            moves = find_move(res, config, sequence_len)
            if not moves: continue
            nconfig = apply_some_move(moves, config)
            nE = energy_function(nconfig, seq)
            if calc_acceptance(cE,nE,seq,0.5):
                cE = nE
                config = nconfig
                if cE<init_energy:
                    #print_at("MC steps beat it...i %d res %d" % (i,res),2)
                    return (False, seq)
    print_at("1k MC steps still failed",3)
    # a little GA
    from minimization import birth, do_mutations, do_cross_over, calc_genetic_energy, choose_fit
    from genetic_lattice import encode_lattice
    indivs = 100
    population = birth(sequence_len, indivs) + [encode_lattice(orig_config)]
    for gen in xrange(50):
        do_cross_over(population,1.0,sequence_len)
        do_mutations(population,1.0,sequence_len)
        population = choose_fit(population, indivs, seq)
        lowest_V = calc_genetic_energy(population[0],seq)
        if lowest_V<init_energy:
            print_at("initial energy %s" % init_energy,4)
            print_at("found lower energy with GA... %s" % lowest_V,5)
            return (False, seq)
    print_at("local minimum? %s" % seq,6)
    print "\n\n\n"
    return (True, seq)


def check_global_minimum(seq, orig_config):
    from genetic_lattice import print_at
    print_at("(global) checking... %s" % seq,7)
    f = open("/tmp/global.out",'a')
    f.write(str(seq))
    f.close()
    # 1mil MC steps
    from minimization import apply_some_move, find_move, calc_acceptance
    from energy_function import energy_function
    config = orig_config
    init_energy = energy_function(config, seq)
    from minimization import calc_acceptance
    sequence_len = len(seq)
    sequence_range = xrange(sequence_len)
    cE = energy_function(config, seq)
    for i in xrange(10**6):
        for res in sequence_range:
            moves = find_move(res, config, sequence_len)
            if not moves: continue
            nconfig = apply_some_move(moves, config)
            nE = energy_function(nconfig, seq)
            if calc_acceptance(cE,nE,seq,0.5):
                cE = nE
                config = nconfig
                if cE<init_energy:
                    print_at("MC steps beat it...i %d res %d" % (i,res),2)
                    return (False, seq)
    print "10k MC steps still failed"
    # a little GA
    from minimization import birth, do_mutations, do_cross_over, calc_genetic_energy, choose_fit
    from genetic_lattice import encode_lattice
    indivs = 100
    population = birth(sequence_len, indivs) + [encode_lattice(orig_config)]
    for gen in xrange(50):
        do_cross_over(population,1.0,sequence_len)
        do_mutations(population,1.0,sequence_len)
        population = choose_fit(population, indivs, seq)
        lowest_V = calc_genetic_energy(population[0],seq)
        if lowest_V<init_energy:
            print "initial energy %s" % init_energy
            print "found lower energy with GA... %s" % lowest_V
            return (False, seq)
    print "global minimum!!!"
    return (True, seq)



# GA stuff

from energy_function import energy_function
from random import randrange, random, choice, sample
randdir = lambda:randrange(6)

def birth(sequence_len, num):
    def live_birth():
        baby = []
        neck,shoulders = -1,-1
        head = -1
        for i in xrange(sequence_len-1):
            while head==neck or head==shoulders:
                head = randdir()
            baby.append(head)

            neck = head
            shoulders = neck ^ 1 # odd<->even
        return baby

    return [live_birth() for j in xrange(num)]


def do_mutations(population, rate, sequence_len):
    i=0
    lenpop = len(population)
    mutated = []
    while i<lenpop:
        i+=1
        if random() > rate: continue

        indiv = list(choice(population))

        mutate_at = randrange(sequence_len-1)
        indiv[mutate_at] = randdir()
        mutated.append(indiv)

    population.extend(mutated)

def do_cross_over(population, rate, sequence_len):
    i=0
    lenpop = len(population)
    crossed = []
    while i<lenpop:
        i+=1
        if random() > rate: continue

        first, second = sample(population,2)
        splice_at_start = randrange(sequence_len-1)
        splice_at_end = randrange(sequence_len-1)

        if splice_at_end<splice_at_start:
            tmp = splice_at_start
            splice_at_start = splice_at_end
            splice_at_end = tmp

        crossed.append(first[:splice_at_start]+second[splice_at_start:splice_at_end]+first[splice_at_end:])
        crossed.append(second[:splice_at_start]+first[splice_at_start:splice_at_end]+second[splice_at_end:])
    population.extend(crossed)

from genetic_lattice import decode_lattice
def calc_genetic_energy(config,sequence):
    return energy_function(decode_lattice(config),sequence)

def choose_fit(population,num,sequence,resort=True):
    if resort:
        population.sort(key=lambda c:calc_genetic_energy(c,sequence))
    if len(population)>num:
        return population[:num]
    return population


# MC stuff

from montecarlo_ss import find_move,calc_acceptance

def apply_move(move, config):
    #move = choice(moves)
    new_config = list(config)
    for i in range(len(move[0])):
        new_config[move[0][i]]=move[1][i]
    return new_config

def apply_some_move(moves, config):
    move = choice(moves)
    new_config = list(config)
    for i in range(len(move[0])):
        new_config[move[0][i]]=move[1][i]
    return new_config
