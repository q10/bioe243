class Point(object):
    def __init__(self, x, y, z):
        self.x=x
        self.y=y
        self.z=z

def energy_function(configuration, sequence):
    i=0
    j=0
    V=0
    lconfig = len(configuration)
    all_neighbors = count_all_neighbors(configuration, lconfig)
    print all_neighbors 
    while i<lconfig:
        pos_i = configuration[i]
        #~ print 'x = ' + str(pos_i[0])
        #~ print 'x = ' + pos_i.x
        bead_type_i = sequence[i]
        neighbors_i = all_neighbors[i]
        lambda_i = calc_lambda(bead_type_i, i, configuration, lconfig, neighbors_i)
        #~ N = count_neighbors(i, configuration, lconfig)
        if all_neighbors[i]!=0:
            j=i+2
            for j in range(j, lconfig-2):
                pos_j = configuration[j]
                bead_type_j = sequence[j]
                neighbors_j = all_neighbors[j]
                deltaX = pos_i[0] - pos_j[0]
                deltaY = pos_i[1] - pos_j[1]
                deltaZ = pos_i[2] - pos_j[2]
                r_squared=((deltaX*deltaX)+(deltaY*deltaY)+(deltaZ*deltaZ))
                #~ print 'r = ' + str(r_squared)
                if r_squared == 1:
                    #~ print 'Found neighbor! Bead i:j = ' + str(i) + ':' + str(j) + '. Positions: ' + str(pos_i) + str(pos_j)
                    lambda_j = calc_lambda(bead_type_j, j, configuration, lconfig, neighbors_j)
                    lambda_avg = 0.5 * (lambda_i + lambda_j)

#~              Interactions between particles
#~              HH -> Bu = Bf = -1
#~              HP -> Bu = 0 and Bf = 1
#~              PH -> Bu = 0 and Bf = 1
#~              PP -> Bu = 1 and Bf = -1
#~              V = (1-l)Bu + l*Bf

                    if (bead_type_i == 'H' and bead_type_j == 'H'):
                        Bu = Bf = -1
                        V = V + ( (1-lambda_avg)*Bu + lambda_avg*Bf )
                    if (bead_type_i == 'H' and bead_type_j == 'P') or (bead_type_i == 'P' and bead_type_j == 'H'):
                        Bu=0
                        Bf=1
                        V = V + ( (1-lambda_avg)*Bu + lambda_avg*Bf )
                    if (bead_type_i == 'P' and bead_type_j == 'P'):
                        Bu=1
                        Bf=-1
                        V = V + ( (1-lambda_avg)*Bu + lambda_avg*Bf )
            j+=1
        i+=1
    return V

def calc_lambda(bead_type, bead_number, configuration, lconfig, neighbors):
    l=0             #~ l = lambda
    if neighbors >= 3:
        l=1
        exit
    si=neighbors
    if bead_type == 'P':
        si0=2
        l=si / si0
    else:
        si0=3
        l=si / si0
    return l

def count_neighbors(bead_number, configuration, lconfig):
    number_of_neighbors=0
    bead_position_i=configuration[bead_number]
    j=0
    while j<lconfig:
        if (abs(bead_number-j)>1):
            bead_position_j = configuration[j]
            deltaX = bead_position_i[0] - bead_position_j[0]
            deltaY = bead_position_i[1] - bead_position_j[1]
            deltaZ = bead_position_i[2] - bead_position_j[2]
            r_squared=((deltaX*deltaX)+(deltaY*deltaY)+(deltaZ*deltaZ))
            if r_squared == 1:
                number_of_neighbors+=1
        j+=1
    #~ print str(bead_number) + ' has ' + str(number_of_neighbors) + ' neighbors.\n'
    return number_of_neighbors

def count_all_neighbors(configuration, lconfig):
    neighbor_space = []
    number_of_neighbors = 0
    i=0
    while i<lconfig:
        number_of_neighbors=0
        bead_position_i = configuration[i]
        j=0
        while j<lconfig:
            if (abs(i-j)>1):
                bead_position_j = configuration[j]
                deltaX = bead_position_i[0] - bead_position_j[0]
                deltaY = bead_position_i[1] - bead_position_j[1]
                deltaZ = bead_position_i[2] - bead_position_j[2]
                r_squared=((deltaX*deltaX)+(deltaY*deltaY)+(deltaZ*deltaZ))
                if r_squared == 1:
                    number_of_neighbors+=1
                    #~ print str(i) + ' neighbors ' + str(j) + '. Number of neighbors = ' +str(number_of_neighbors) + '\n'
            j+=1
        neighbor_space.append(number_of_neighbors)
        #~ print str(i) + ' has ' + str(number_of_neighbors) + ' neighbors.\n'
        i+=1
    #~ print str(bead_number) + ' has ' + str(number_of_neighbors) + ' neighbors.\n'
    #~ print neighbor_space
    return neighbor_space

if __name__ == '__main__':
    filename = 'C:\Users\Kristin\Downloads\protein.36'
    #{
    import sys
    if len(sys.argv)>1:
        filename = sys.argv[1]
    #}
    configuration = []
    r_file = open(filename, 'r')
    r_line = r_file.readline()
    while r_line != '':
        r = map(float, r_line.strip().replace(',', ' ').split())
        #r = Point(coord[0],coord[1],coord[2])
        configuration.append(r)
        #~ print configuration
        r_line = r_file.readline()
    r_file.close()
    sequence = 'HHPHHPPHPPPHPPHPHHHHPPPPHHPPHPHHHHPP'
    V = energy_function(configuration,sequence)
    print 'Seq 1: V = ' +str(V)
    sequence = 'HHHPPHHPHHHHPHHPPHPHHPHPPHPPPHHPPPPP'
    V = energy_function(configuration,sequence)
    print 'Seq 2: V = ' +str(V)


#~ Conventions to use:
#~ Underscores rather than camel case
