#####################
## Energy Function ##
#####################


# long explanation for these weird numbers, but we can't use Point for them
nbr_gen_3d = lambda i,p,nbr: p + ((1-((nbr&1)<<1)) if (((nbr&6)>>1)==i) else 0)

#get_neighbors = lambda pos: [tuple(nbr_gen_3d(i,pos[i],nbr) for i in (0,1,2)) for nbr in (0,1,2,3,4,5)]
# memoized version of get_neighbors since it gets called so often
neighbor_memoizations = {}
# assumes no lattice goes more than 35 to one direction (can't for us)
nbr_hash = lambda pos: pos[0]+35*pos[1]+1225*pos[2]
def get_neighbors(pos):
    hash = nbr_hash(pos)
    if hash not in neighbor_memoizations:
        nbrs = [tuple(nbr_gen_3d(i,pos[i],nbr) for i in (0,1,2)) for nbr in (0,1,2,3,4,5)]
        neighbor_memoizations[hash] = nbrs
        return nbrs
    return neighbor_memoizations[hash]

get_full = lambda nbrs,config: [nbr for nbr in nbrs if nbr in config]
get_full_neighbors = lambda i,config: get_full(get_neighbors(config[i]),config)

# H=0, P=1
convert_sequence = lambda sqn: tuple(0 if r=='H' else 1 for r in sqn)
#Bu_model = ((-1,0),(0,1))
Bu_model = lambda i,j: i+j-1
#Bf = ((-1,1),(1,-1))
# optimization
Bf_model = lambda i,j: -1 if i==j else 1

seen_lambdas = []
def get_lambda_i(i,residue_i,neighbor_len, sequence_len):
    seen_len = len(seen_lambdas)
    if seen_len>i: return seen_lambdas[i]
    s_i = neighbor_len - 2
    if i==0 or i==sequence_len-1:
        s_i += 1
    s_i0 = 3 - residue_i
    lambda_i = float(s_i)/s_i0
    if lambda_i>1: lambda_i=1.0

    # check prolly isn't necessary
    if seen_len==i: seen_lambdas.append(lambda_i)
    return lambda_i

def has_repeats(config):
    seen = {}
    for item in config:
        #print item
        if item in seen: return True
        seen[item] = 1
    return False

def energy_function(configuration, sequence):
    if has_repeats(configuration):
        #print 'has repeats... %s' % configuration
        return 0
    #print 'no repeats...'
    del seen_lambdas[:] # clear it

    i=0
    V_lattice = 0
    lconfig = len(configuration)
    neighbor_space = [get_full_neighbors(j,configuration) for j in xrange(lconfig)]
    neighbor_len_space = map(len, neighbor_space)
    #print neighbor_len_space
    #print neighbor_space
    #print len(neighbor_space)
    sequence_len = len(sequence)
    while i < lconfig:
        residue_i = sequence[i]
        nbrs = neighbor_space[i]
        lambda_i = get_lambda_i(i,residue_i,neighbor_len_space[i],sequence_len)
        for nbr in nbrs:
            nbrj = configuration.index(nbr)
            if nbrj<=i: continue

            # skip bonded neighbors (already skipped if using non_bonded_neighbors)
            if nbrj-i==1 or nbrj-i==-1: continue

            residue_j = sequence[nbrj]
            Bu = Bu_model(residue_i,residue_j)
            Bf = Bf_model(residue_i,residue_j)
            lambda_j = get_lambda_i(nbrj,residue_j,neighbor_len_space[nbrj],sequence_len)

            lambdaij = (lambda_i + lambda_j)/2.0
            Vij = (1-lambdaij)*Bu+lambdaij*Bf
            V_lattice += Vij
        i+=1
    return V_lattice

if __name__=='__main__':
    filename = 'protein.36'
    #{
    import sys
    if len(sys.argv)>1:
        filename = sys.argv[1]
    #}
    configuration = []
    lines = open(filename, 'r').readlines()
    for line in lines:
        pos = line.strip().replace(',',' ').split(' ')
        configuration.append(tuple(int(p) for p in pos))
    sequence = convert_sequence('HHPHHPPHPPPHPPHPHHHHPPPPHHPPHPHHHHPP')
    V = energy_function(configuration,sequence)
    print 'Seq 1: V = ' +str(V)
    sequence = convert_sequence('HHHPPHHPHHHHPHHPPHPHHPHPPHPPPHHPPPPP')
    V = energy_function(configuration,sequence)
    print 'Seq 2: V = ' +str(V)


