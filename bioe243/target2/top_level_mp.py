

# need 5 global minima...


# ::::::::::::::::::::::::::WINS:::::::::::::::::::::::::::
# multiprocessing way:
# start with an empty list of results
# loop cpu_count() times
#   create a random sequence
#   create a apply_async job to check local minimum, store asyncresults
# infinite loop
#   loop through local asyncresults
#       if ready, get result
#           if result means local minimum
#               creaty apply_async job to check global minimum and store seperately
#   loop through global asyncresults ( replaced by callback from local )
#       if ready, get result
#           if result means global minimum
#               store to list of results
#               if we have 5 global minimum, tell pool to terminate and break
# we have 5 global minimum if we reached here

from multiprocessing import Pool,cpu_count
from random_seq import fastest_create_random_sequence as create_random_sequence
from minimization import find_move
from minimization import check_global_minimum, check_local_minimum

def calc_moves(config):
    sequence_len = len(config)
    moves = []
    for i in range(sequence_len):
        moves.extend(find_move(i,config,sequence_len))
    return moves

to_try = [
    [0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0]
]

class GlobalChecker(object):
    def __init__(self, pool, num_to_find, config):
        self.global_minima = []
        self.pool = pool
        self.local_results = []
        self.still_checking = True
        self.num_to_find = num_to_find
        self.config = config
        self.moves = calc_moves(config)

    def start(self, result):
        if result[0]: # local minimum :)
            seq = result[1]
            self.pool.apply_async(check_global_minimum, (seq,self.config), callback=self.done)

    def done(self, result):
        if result[0]: # global minimum :)
            seq = result[1]
            self.global_minima.append(seq)
            print "11 %d %d" % (len(self.global_minima),self.num_to_find)
            if len(self.global_minima)>=self.num_to_find:
                print "8"
                self.output_sequences_to_file()
                print "9"
                self.still_checking = False
                self.pool.join()#.terminate()
                print "10"

    def output_sequences_to_file(self):
        f = open('/tmp/global_minimum.output','a')
        f.write('\n\n\n%s\n\n\n' % str(self.global_minima))
        f.close()

    def check_local(self):
        #if len(to_try)>0:
        #    seq = to_try.pop()
        #    print "trying %s" % seq
        #else:
        #    seq = create_random_sequence()
        seq = create_random_sequence()
        if self.still_checking:
            result = self.pool.apply_async(check_local_minimum, (seq,self.moves,self.config), callback=self.start)
            self.local_results.append(result)

    def restart_locals(self):
        #print "5 %d" % len(self.local_results)
        i = 0
        curlen = len(self.local_results)
        while i < curlen:
            #print "6 %d" % i
            result = self.local_results[i]
            if not result.ready():
                result.wait(1)
                i += 1
            else:
                self.local_results.pop(i)
                curlen -= 1
                # post-processing handled by callback
                # restart process if we're not done
                if self.still_checking:
                    self.check_local()


def main(num_to_find=5):
    from genetic_lattice import decode_lattice, decode_string
    config = decode_lattice(decode_string('RBBBRFFFUBUFRDDBUBLLLDDFUUFRB'))
    pool = Pool()
    global_checker = GlobalChecker(pool, num_to_find, config)
    print "1"
    for i in range(cpu_count()*30):
        global_checker.check_local()
    print "3"

    #i = 0
    while True:
        #print "4 %d" % i
        #i += 1
        #if i>10: break
        global_checker.restart_locals()

    #pool.join()
    from time import sleep
    sleep(5)

# possible ones from previous runs:
# [0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1]



if __name__ == '__main__':
    #from random_seq import seen_seq
    #f = open('/tmp/seen_seq.out','r')
    #for line in f.xreadlines():
    #    if line=='': break
    #    seen_seq.append(int(line))
    #f.close()
    try:
        main(1)
    finally:
        #from random_seq import seen_seq
        #f = open('/tmp/seen_seq.out','w')
        #for s in seen_seq:
        #    f.write(str(s)+'\n')
        #f.close()
        pass



# rejected ways:

# normal way:
# start with an empty list of them
# infinite loop
#   create a random sequence
#   check if it's a local minimum
#   if not, continue
#   else, check if it's a global minimum
#       if global minimum
#           store in list
#           if we have 5 global minimum, exit
# we have 5 global minimum if we reached here

# pp way: too complicated...
# start with an empty list of them
# create job template for local and global minimum checking
#   set callbacks
# loop get_ncpus() times
#   create a random sequence
#   submit template for sequence
# callback for local minimum checking
#   if its a local minimum, submit template for global
#   else start a blah

