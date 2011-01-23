# python global_minimum.py [first|second] [# indivs per gen] [# generations]

# a huge speedup for GA would be to use python-blist or another copy-on-write list method

from energy_function import convert_sequence, energy_function
from genetic_lattice import decode_lattice,encode_lattice
from heat_capacity import MC_step

from random import randrange, random, sample, choice
import sys

randdir = lambda:randrange(6)

def birth(sequence_len, num):
    #return [[randdir() for i in xrange(sequence_len-1)] for j in xrange(num)]
    # try to create a more likely valid config by biasing against easily-invalid ones
    # for instance, any config with [0,1] is invalid because it overlaps right away
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

def do_cross_over(population, rate, sequence_len):
    #while random() < rate:
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

def calc_genetic_energy(config,sequence):
    return energy_function(decode_lattice(config),sequence)

# two ways to choose fit:
#   shedskin a slow, non-closure-based version
#   or use a pythonic, closure-based version
def choose_fit(population,num,sequence,resort=True):
    if resort:
        population.sort(key=lambda c:calc_genetic_energy(c,sequence))
    #population = sorted(population,key=lambda c:calc_genetic_energy(c,sequence))
    #custom_selection_sort(population,sequence) # for shedskin
    if len(population)>num:
        return population[:num]
    return population

print_cursor = 0
import os
ttycolumns = int(os.popen('stty size', 'r').read().split()[1])
def print_at(st,line=0):
    global print_cursor
    st = str(st)
    if line>print_cursor:
        sys.stdout.write('\n' * (line-print_cursor))
    elif line<print_cursor:
        sys.stdout.write('\033[%dA' % (print_cursor-line))

    num_lines = (len(st)/ttycolumns)
    print_cursor = line+1+num_lines
    #print len(st),ttycolumns,st
    print st+' '*(ttycolumns-len(st)%ttycolumns)

from time import time
class TimingContext(object):
    def __init__(self,at):
        self.line = at
    def __enter__(self):
        self.start = time()
    def __exit__(self,a,b,c):
        print_at("Took %f" % (time()-self.start),self.line)

def global_minimum():
    if len(sys.argv)==1 or 'first' in sys.argv:
        sequence = convert_sequence('HHPHHPPHPPPHPPHPHHHHPPPPHHPPHPHHHHPP')
    else:
        sequence = convert_sequence('HHHPPHHPHHHHPHHPPHPHHPHPPHPPPHHPPPPP')
    sequence_len = len(sequence)

    if len(sys.argv)>1 and not sys.argv[1].isdigit(): sys.argv.pop(1) # remove first/second

    lenargv = len(sys.argv)
    # GA data:
    # number of individuals per generation
    ind_per_gen  = int(sys.argv[1]) if lenargv>1 else 100
    # number of generations
    num_generations = int(sys.argv[2]) if lenargv>2 else 100

    # number of GA/SA swaps
    num_tries = int(sys.argv[3]) if lenargv>3 else 10

    # SA data:
    alpha = 0.995
    ntemp = int(sys.argv[4]) if lenargv>4 else 1000
    ncycles = int(sys.argv[5]) if lenargv>5 else 100
    orig_T_star = float(sys.argv[6]) if lenargv>6 else 7.0
    # T_star should be <0.2 after ntemp, so
    #    T_star*alpha**ntemp<0.2
    #    (0.2/T_star)**(1.0/ntemp)=alpha
    alpha = (0.2/orig_T_star)**(1.0/(ntemp+1))
    print alpha


    # mutation rates
    rate_cross = 0.8
    rate_mutate = 0.8

    population = birth(sequence_len,ind_per_gen)

    lowest_V = lowest_V_sofar = 0
    highest_V = 1
    equal_generations = 0
    lower_generations = []
    for numtry in xrange(num_tries):
        choose_fit(population,ind_per_gen,sequence) # just so it sorts
        for i in xrange(num_generations):
            print_at('generation %d' % i,0)
            # mate best half?
            with TimingContext(1):
                #mating_population = choose_fit(population,ind_per_gen/2,sequence,resort=False)
                mating_population = population[:ind_per_gen/2]
                population = population[ind_per_gen/2:] # population is pre-sorted by choose_fit

            with TimingContext(2):
                if lowest_V < -1.0 and lowest_V == highest_V: # rebirth soon necessary
                    equal_generations += 1
                else:
                    equal_generations = 0

                if equal_generations > 4: # we lost genetic diversity...
                    mating_population = birth(sequence_len,ind_per_gen)
                    population = sample(population,2)


            with TimingContext(3):
                do_cross_over(mating_population,rate_cross,sequence_len)

            with TimingContext(4):
                do_mutations(mating_population,rate_mutate,sequence_len)




            # keep best ones only
            with TimingContext(5):
                population.extend(mating_population)
            with TimingContext(6):
                population = choose_fit(population,ind_per_gen,sequence)

            # select lowest energy
            with TimingContext(7):
            #if True:
                lowest_conf = population[0]
                lowest_V = calc_genetic_energy(lowest_conf,sequence)
                highest_V = calc_genetic_energy(population[-1],sequence)

            if lowest_V<lowest_V_sofar:
                lowest_V_sofar = lowest_V
                print_at("lowest V: %f for config %s" % (lowest_V,str(lowest_conf).replace(' ','')),8)
                lower_generations.append((i,lowest_V,lowest_conf))
            print_at("current_lowest %f" % lowest_V,9)
            print_at("current_highest %f" % highest_V,10)
            print_at("Generations worth their CPU: %s" % str(lower_generations).replace(' ',''), 11)

        # have lowest conf, try annealing now...
        lowest_conf = decode_lattice(lowest_conf)
        T_star = orig_T_star
        for k in xrange(ntemp):
            print_at("SA temperature %f (%d)" % (T_star,k),25)
            current_config = lowest_conf
            current_V = lowest_V
            print_at("lowest_V: %f" % lowest_V,26)
            for j in xrange(ncycles):
                print_at("Cycle %d" % j,27)
                next_conf,next_V = MC_step(current_config, sequence, T_star)
                if next_V<lowest_V:
                    current_V = next_V
                    current_config = next_conf
                print_at("current V: %f" % next_V,28)
            # take best config from cycles and use that
            if current_V < lowest_V:
                lowest_V = current_V
                lowest_conf = current_config
                print_at("lowest V: %f for config %s" % (lowest_V,str(encode_lattice(lowest_conf)).replace(' ','')),29)
            T_star *= alpha

        # aaaaand start over
        population = [encode_lattice(lowest_conf)]+birth(sequence_len,ind_per_gen)
        mating_population = []

    print_at('',30)


if __name__ == '__main__':
    #{
    #import cProfile
    #cProfile.run('global_minimum()')
    #sys.exit(0)
    #}
    try:
        global_minimum()
    except:# KeyboardInterrupt:
        import traceback
        print_at(traceback.format_exc(),30)
