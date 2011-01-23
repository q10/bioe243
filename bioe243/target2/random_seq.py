import random
from bisect import bisect_left

#{
def create_random_sequence():
    while True:
        ran = bin(random.getrandbits(30))[2:]
        if ran.count('0') == ran.count('1'):
            ran = ran.replace('0','\0')
            ran = ran.replace('1','\1')
            break
    return ran


range_30 = range(30)
def faster_create_random_sequence():
    # wow, randbits is way faster even though equality isn't guaranteed
    Hpos = random.sample(range_30,15)
    return [0 if i in Hpos else 1 for i in range_30]
#}

orig_seq = [0]*15+[1]*15
def get_seq_int(seq):
    seq_int = 0
    for i in xrange(30):
        if orig_seq[i]:
            seq_int += 1<<i
    return seq_int

#from sortedlist import SortedCollection
#seen_seq = SortedCollection()
seen_seq = []

def check_in(i,key):
    if i != len(seen_seq) and seen_seq[i] == key:
        return True
    return False


def fastest_create_random_sequence_and_check():
    while True:
        random.shuffle(orig_seq)
        seq_int = get_seq_int(orig_seq)
        i = bisect_left(seen_seq, seq_int)
        if not check_in(i,seq_int):
            seen_seq.insert(i,seq_int)
            return orig_seq

def fastest_create_random_sequence():
    random.shuffle(orig_seq)
    return orig_seq

if __name__ == '__main__':
    num = 5000
    #{
    #from time import time
    #start = time()
    for i in xrange(num):
        [1 if ord(i)==1 else 0 for i in create_random_sequence()]
    #print "original took %f" % (time()-start)

    #start = time()
    for i in xrange(num):
        faster_create_random_sequence()
    #print "faster took %f" % (time()-start)

    for i in xrange(num):
        fastest_create_random_sequence_and_check()
    #}

    #start = time()
    for i in xrange(num):
        fastest_create_random_sequence()
    #print "fastest took %f" % (time()-start)
