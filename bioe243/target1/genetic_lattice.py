LEFT,RIGHT,DOWN,UP,BACK,FRONT = 0,1,2,3,4,5
DIRECTIONS = (LEFT,RIGHT,DOWN,UP,BACK,FRONT)

def dir_dim_change(curr,next):
    for i in range(len(curr)):
        if curr[i]!=next[i]:
            dim = i
            if next[i]>curr[i]:
                dir = 1
            else:
                dir = -1
    return dir,dim

def encode_lattice(config):
    encoding = []
    for i in xrange(len(config)-1):
        dir,dim = dir_dim_change(config[i],config[i+1])
        if dir<0: dir = 0
        encoding.append((dim*2)+dir) # automatic
    return encoding

def decode_lattice(encoding):
    config = [(0,0,0)]*(len(encoding)+1)
    #for i in xrange(1,len(encoding)):
    i = 0
    encoding_len = len(encoding)
    while i<encoding_len:
        e = encoding[i]
        dim = e>>1
        dir = ((e&1)<<1)-1
        i += 1
        next_config = tuple([p+dir if d==dim else p for d,p in enumerate(config[i-1])])
        config[i] = next_config
    return config

def choose_fit(population,num):
    population.sort(key=calc_genetic_energy)
    if len(population)>numpop:
        return population[0:numpop]
    return population
