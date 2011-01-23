from genetic_lattice import decode_string, decode_lattice, print_at
from energy_function import energy_function

def sequence_from_i(i,sequence_len):
    return tuple((i>>j)&1 for j in range(sequence_len))

def num_H(sequence):
    return sum(sequence)

if __name__=='__main__':
    configuration = decode_lattice(decode_string('RBBBRFFFUBUFRDDBUBLLLDDFUUFRB'))
    sequence_len = len(configuration)
    sequence_len2 = len(configuration)/2

    max_configs = int(2**sequence_len)
    #start = 10000
    start = 32767 # first half-H config

    #found = []
    lowest_seq = None
    lowest_V = 0
    j = 0
    for i in xrange(start,max_configs):
        sequence = sequence_from_i(i,sequence_len)
        if num_H(sequence)!=sequence_len2: continue
        j += 1
        #found.append(sequence)
        V = energy_function(configuration, sequence)
        if j%10000==0:
            print_at(str(sequence), 2)
            print_at("current V: %f, lowest_V: %f" % (V,lowest_V), 3)
        if V<lowest_V:
            lowest_V = V
            lowest_seq = sequence
            lowest_seq_str = ''.join(['H' if s==0 else 'P' for s in lowest_seq])
            print_at("lowest V: %f for sequence %s" % (lowest_V,lowest_seq_str),4)
    lowest_seq_str = ''.join(['H' if s==0 else 'P' for s in lowest_seq])
    print_at("V: %f for sequence %s" % (lowest_V,lowest_seq_str),5)
    #print_at(len(found),5)
    #print_at(found,6)

