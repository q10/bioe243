def convert_moves(configuration):
    lconfig = len(configuration)
    x = y = z = -1
    for i in range(0,lconfig-1):
        curr = configuration[i]
        next = configuration[i+1]
        for j in range(3):            
            if i == 0 and x == -1:
                dx = next[j] - curr[j]
                if dx != 0:
                    x = j   #define x-axis as the coordinate of the first movement
            if y == -1:
                dy = next[j] - curr[j]
                if dy != 0 and j != x:
                    y = j   #define y-axis as the coordinate for the first non-x movement 	          
            if z == -1:
                dz = next[j] - curr[j]
                if dz != 0 and j != x and j != y:
                    z = j   #define z-axis as the coordinate for the first non-xy movement
    #print x, y, z, dx, dy, dz
    moves = ''
    for k in range(0,lconfig-1):
        curr = configuration[k]
        next = configuration[k+1]
        diff = [m-n for m,n in zip(next,curr)]  #tuple of (next point - current point)
        if diff[x] == dx:
            moves += 'R'   #first movement is in 'R' direction
        elif diff[x] == -dx:
            moves += 'L'
        elif diff[y] == dy:
            moves += 'B'   #first movement out of 1D is in 'B' direction
        elif diff[y] == -dy:
            moves += 'F'
        elif diff[z] == dz:
            moves += 'U'   #first movement off plane is in 'U' direction
        elif diff[z] == -dz:
            moves += 'D'
    print moves     #converted sequence of movements compliant with instructions

if __name__ == '__main__':
    configuration = [(-1,2,0), (-1,3,0), (-1,4,0), (0,4,0), (1,4,0), (1,4,1)]
    convert_moves(configuration)


