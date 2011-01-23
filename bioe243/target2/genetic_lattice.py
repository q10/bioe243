import sys,os
LEFT,RIGHT,BACK,FRONT,DOWN,UP = 0,1,2,3,4,5
DIRECTIONS = (LEFT,RIGHT,BACK,FRONT,DOWN,UP)
TEXTDIR = {'L':0,'R':1,'B':2,'F':3,'D':4,'U':5}

def decode_string(str_enc):
    return [TEXTDIR[s] for s in str_enc]

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

print_cursor = 0
try:
    ttycolumns = int(os.popen('stty size', 'r').read().split()[1])
except:
    ttycolumns = 80
def print_at(thing,line=0):
    global print_cursor
    st = str(thing)
    if line>print_cursor:
        sys.stdout.write('\n' * (line-print_cursor))
    elif line<print_cursor:
        sys.stdout.write('\033[%dA' % (print_cursor-line))

    num_lines = (len(st)/ttycolumns)
    print_cursor = line+1+num_lines
    print st+' '*(ttycolumns-(len(st)%ttycolumns))


