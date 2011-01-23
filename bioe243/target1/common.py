from math import sqrt, log, exp

IMAX = 2.0**31.0

# This code follows exactly the algorithm described in Numerical Recipes
class Ran3:
    m = IMAX-1
    mseed = 161803398

    def __init__(self, init):
        mj = abs(Ran3.mseed-abs(init)) % Ran3.m
        self.rand_table = [0]*55
        self.rand_table[54] = mj
        self.iseed, ii, mk = 0, 0, 1

        # initialize rand_table using init as seed and large number mseed
        for i in range(0, 54):
            ii = (21*i) % 55
            self.rand_table[ii] = mk
            mk = mj - mk
            if mk < 0:
                mk += Ran3.m
            mj = self.rand_table[ii]
        
        # randomize table warming up the generator
        for k in range(0, 4):
            for i in range(0, 55):
                self.rand_table[i] -= self.rand_table[1+((i+30)%54)]
                if (self.rand_table[i] < 0):
                    self.rand_table[i] += Ran3.m

        

    # generate the next random number
    def generate(self):
        self.rand_table[self.iseed] = (self.rand_table[self.iseed] - self.rand_table[(self.iseed+31) % 55])
        if self.rand_table[self.iseed] < 0:
            self.rand_table[self.iseed] += Ran3.m
        self.iseed = (self.iseed+1) % 55

        return self.rand_table[self.iseed] / Ran3.m


class BoxMueller:
    def __init__(self, seed):
        self.ran3 = Ran3(seed)

    def generate(self):
        r = 1.0
        while r >= 1.0:
            x1 = (2.0*self.ran3.generate()) - 1.0
            x2 = (2.0*self.ran3.generate()) - 1.0
            r = (x1**2) + (x2**2)

        r = sqrt(-2.0*log(r)/r)
        return (x1*r, x2*r)



class Bead:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z
        self.energy = 0

    def __str__(self):
        position = "Position " + str(self.position()) + "\n"
        return position + "Energy = " + str(self.energy)
    def position(self):
        return (self.x, self.y, self.z)
    def add_position(self, dx, dy, dz):
        self.x += dx
        self.y += dy
        self.z += dz
    def set_position(self, x, y, z):
        (self.x, self.y, self.z) = (x, y, z)
    def add_energy(self, dE):
        self.energy += dE
    def set_energy(self, e):
        self.energy = e

class AABead(Bead):
    def __init__(self, x, y, z, typ):
        Bead.__init__(self, x, y, z)
        self.type = typ

class Particle(Bead):
    def __init__(self, x, y, z):
        Bead.__init__(self, x, y, z)
        self.vx, self.vy, self.vz = 0, 0, 0
        self.fx, self.fy, self.fz, = 0, 0, 0
        self.boundary = 0
        self.m = 1
        self.is_closed_interval = False
        
    def __str__(self):
        position = "Position " + str(self.position()) + "\n"
        velocity = "Velocity " + str(self.velocity()) + "\n"
        force = "Force "+ str(self.force()) + "\n"
        return position + velocity + force + "Energy = " + str(self.energy)
    
    def velocity(self):
        return (self.vx, self.vy, self.vz)
    def force(self):
        return (self.fx, self.fy, self.fz)

    def add_position(self, dx, dy, dz):
        Bead.add_position(self, dx, dy, dz)
        if self.is_closed_interval:
            self._enclose_to_boundary()
    def add_velocity(self, dvx, dvy, dvz):
        self.vx += dvx
        self.vy += dvy
        self.vz += dvz
    def add_force(self, dfx, dfy, dfz):
        self.fx += dfx
        self.fy += dfy
        self.fz += dfz

    def set_position(self, x, y, z):
        Bead.set_position(self, x, y, z)
        if self.is_closed_interval:
            self._enclose_to_boundary()
    def set_velocity(self, vx, vy, vz):
        (self.vx, self.vy, self.vz) = (vx, vy, vz)
    def set_force(self, fx, fy, fz):
        (self.fx, self.fy, self.fz) = (fx, fy, fz)

    def set_boundary(self, b):
        self.boundary = b
        self.is_closed_interval = True
    
    def _enclose_to_boundary(self):
        if self.x < -self.boundary:
            self.x %= self.boundary
        if self.x > self.boundary:
            self.x %= -self.boundary
        if self.y < -self.boundary:
            self.y %= self.boundary
        if self.y > self.boundary:
            self.y %= -self.boundary
        if self.z < -self.boundary:
            self.z %= self.boundary
        if self.z > self.boundary:
            self.z %= -self.boundary



class System:
    def __init__(self, filename, temp, density):
        self.particles = []
        self.__create_particles(filename)

        # set system parameters
        self.temperature = temp
        self.N = len(self.particles)
        #self.rho = density
        #self.cube_boundary = ((self.N/self.rho)**(1.0/3.0))/2.0
        
        #for particle in self.particles:
        #    particle.set_boundary(self.cube_boundary)

        self._update_energies()
        self._update_forces()
        
    def __create_particles(self, filename):
        # initialize particle positions
        r_file = open(filename, 'r')
        r_line = r_file.readline()
        while r_line != '':
            r = map(float, r_line.strip().replace(',', ' ').split())
            self.particles.append(Particle(r[0],r[1],r[2]))
            r_line = r_file.readline()
        r_file.close()

    def collect_positions(self):
        positions = []
        for particle in self.particles:
            positions.append(particle.position())
        return positions

    def collect_velocities(self):
        velocities = []
        for particle in self.particles:
            velocities.append(particle.velocity())
        return velocities

    def collect_forces(self):
        forces = []
        for particle in self.particles:
            forces.append(particle.force())
        return forces

    def collect_energies(self):
        energies = []
        for particle in self.particles:
            energies.append(particle.energy)
        return energies

    def _update_forces(self):
        for particle in self.particles:
            particle.set_force(0,0,0)
        for i in range(self.N-1):
            for j in range(i+1, self.N):
                (dfx, dfy, dfz) = self.force_between(self.particles[i], self.particles[j])    
                self.particles[i].add_force(-dfx, -dfy, -dfz)
                self.particles[j].add_force(dfx, dfy, dfz)

    def force_between(self, first, other):
        dx, dy, dz = first.x-other.x, first.y-other.y, first.z-other.z
        #(dx, dy, dz) = self.__mirror_convention(first.x-other.x, first.y-other.y, first.z-other.z)
        r = sqrt((dx**2.0)+(dy**2.0)+(dz**2.0))
        tmp_force = -24.0 * ((2.0 * ((1/r)**14.0)) - ((1/r)**8.0))
        return (tmp_force*dx, tmp_force*dy, tmp_force*dz)

    def _update_energies(self):
        for particle in self.particles:
            particle.set_energy(0)
        for i in range(self.N-1):
            for j in range(i+1, self.N):
                dE = self.energy_between(self.particles[i], self.particles[j])
                self.particles[i].add_energy(dE)
                self.particles[j].add_energy(dE)

    def energy_between(self, first, other):
        dx, dy, dz = first.x-other.x, first.y-other.y, first.z-other.z
        #(dx, dy, dz) = self.__mirror_convention(first.x-other.x, first.y-other.y, first.z-other.z)
        r = sqrt((dx**2.0)+(dy**2.0)+(dz**2.0))
        return 4.0 * (((1/r)**12.0) - ((1/r)**6.0))

    def _restore_old_energies(self, old_energies_table):
        for i in range(self.N):
            self.particles[i].set_energy(old_energies_table[i])

    def _restore_old_positions(self, old_positions_table):
        for i in range(self.N):
            (x, y, z) = old_positions_table[i]
            self.particles[i].set_position(x, y, z)

    def _restore_old_velocities(self, old_velocities_table):
        for i in range(self.N):
            (vx, vy, vz) = old_velocities_table[i]
            self.particles[i].set_velocity(vx, vy, vz)

    def _restore_old_forces(self, old_forces_table):
        for i in range(self.N):
            (fx, fy, fz) = old_forces_table[i]
            self.particles[i].set_force(fx, fy, fz)

    # Enforce mirror image convention
    def __mirror_convention(self, dx, dy, dz):
        if dx < -self.cube_boundary:
            dx %= self.cube_boundary
        if dx > self.cube_boundary:
            dx %= -self.cube_boundary
        if dy < -self.cube_boundary:
            dy %= self.cube_boundary
        if dy > self.cube_boundary:
            dy %= -self.cube_boundary
        if dz < -self.cube_boundary:
            dz %= self.cube_boundary
        if dz > self.cube_boundary:
            dz %= -self.cube_boundary
        return (dx, dy, dz)

    def system_potential_energy(self):
        return sum(self.collect_energies())/2.0
    def system_forces(self):
        return self.collect_forces()
    def system_positions(self):
        return self.collect_positions()



class Lattice:
    def __init__(self, filename, temp):
        self.beads = []
        self.__create_beads(filename)

        # set system parameters
        self.temperature = temp
        self.N = len(self.beads)
        self.rand = Ran3(234257424417)
        
        self._update_energies()
        
    def __create_beads(self, filename):
        # initialize bead positions
        r_file = open(filename, 'r')
        r_line = r_file.readline()
        while r_line != '':
            r = r_line.strip().replace(',', ' ').split()
            #~ typ = r.pop(0)
            r = map(float, r)
            self.beads.append(AABead(r[0], r[1], r[2]))
            r_line = r_file.readline()
        r_file.close()

    def _move_bead_to(self, i, coord):
        (x, y, z) = coord
        self.beads[i].set_position(x, y, z)

    def _update_energies(self):
        for bead in self.beads:
            bead.set_energy(self.energy_of_bead(bead))
    
    def energy_of_bead(self, bead):
        energy = 0.0
        neighbor_coords = self.neighbor_coords_of_bead(bead)
        for coord in neighbor_coords:
            possible_bead = self.find_bead_at_this_coordinate(coord)
            energy += self.energy_between(bead, possible_bead)
        return energy

    def collect_positions(self):
        positions = []
        for bead in self.beads:
            positions.append(bead.position())
        return positions
    
    def collect_energies(self):
        energies = []
        for bead in self.beads:
            energies.append(bead.energy)
        return energies

    def system_potential_energy(self):
        return sum(self.collect_energies())/2.0
    
    def _restore_old_positions(self, old_positions_table):
        for i in range(self.N):
            (x, y, z) = old_positions_table[i]
            self.beads[i].set_position(x, y, z)

    def _restore_old_energies(self, old_energies_table):
        for i in range(self.N):
            self.beads[i].set_energy(old_energies_table[i])

    def save_old_system(self):
        return(self.collect_positions(), self.collect_energies())

    def _restore_old_system(self, old_system_table):
        self._restore_old_positions(old_system_table[0])
        self._restore_old_energies(old_system_table[1])

    def neighbor_coords_of_bead(self, bead):
        (x, y, z) = bead.position()
        neighbor_coords = []
        neighbor_coords.append((x+1,y,z))
        neighbor_coords.append((x-1,y,z))
        neighbor_coords.append((x,y+1,z))
        neighbor_coords.append((x,y-1,z))
        neighbor_coords.append((x,y,z+1))
        neighbor_coords.append((x,y,z-1))
        return neighbor_coords

    def find_bead_at_this_coordinate(self, coords):
        for bead in self.beads:
            if bead.position() == coords:
                return bead
        return None
    
    def energy_between(self, first, other, energy_model=1):
        energy = 0.0;
        # if using energy model 0 (w/o solvent) and index is legit
        if energy_model == 0 and other != None:
            if first.type == 'H':
                if other.type == 'H':
                    energy -= 1.0
            else:
                if other.type == 'P':
                    energy -= 1.0
          
        # else we are using energy model 1 (w/ solvent)
        else:
            if other != None:
                if first.type == 'H':
                    if other.type == 'H':
                        energy -= 1.0
                else:
                    if other.type != 'H':
                        energy -= 0.5
            else:
                if first.type == 'H':
                    energy += 0.5;
                else:
                    energy -= 0.5;

        return energy

    def corner_flip_possible_for_bead(self, i):
        # account for end molecules
        if i == 0 or i == self.N-1:
            return []
        possible_corner = self.corner_coord_of_bead(i)
        if possible_corner == None:
            return []
        elif self.find_bead_at_this_coordinate(possible_corner) == None:
            return [possible_corner]
        return []

    def corner_coord_of_bead(self, i):
        # account for end molecules
        if i == 0 or i == self.N-1:
            return None
        
        left = self.beads[i-1]
        middle = self.beads[i]
        right = self.beads[i+1]
        
        # check if left, middle, and center are NOT colinear
        code = 0
        if middle.x == left.x and middle.x == right.x:
            code += 1
        if middle.y == left.y and middle.y == right.y:
            code += 1
        if middle.z == left.z and middle.z == right.z:
            code += 1
        if code > 1:
            return None
                
        # find the opposite corner's coordinates
        corner = None
        if middle.x == left.x and middle.x == right.x:
            if middle.y == left.y:
                corner = (middle.x, right.y, left.z)
            else:
                corner = (middle.x, left.y, right.z)
        elif middle.y == left.y and middle.y == right.y:
            if middle.x == left.x:
                corner = (right.x, middle.y, left.z)
            else:
                corner = (left.x, middle.y, right.z)
        elif middle.z == left.z and middle.z == right.z:
            if middle.x == left.x:
                corner = (right.x, left.y, middle.z)
            else:
                corner = (left.x, right.y, middle.z)
        else:
            return None
        
        return corner
            
    def crankshafts_possible_for_bead(self, index):
        '''
            The crankshaft will be referred to as follows
                             (c = index)____________(d)
                                |                    |
                                |                    |
                                |                    |
                      (a)______(b)                  (e)______(f)

            returns the number of crankshafts possible with amino acid index and index+1 
            (just the ones that are both open and possible to rotate to), and 
            places the new coords in pairs into available_next_coords' row 0-5      
        '''
        # crankshaft cannot be performed on the first or last two amino acids
        if index < 2 or index > self.N-3:
            return []

        # if this position is not a crankshaft to begin with, then return null
        all_crankshafts_possible = self.all_crankshafts_possible_for_bead(index)
        if all_crankshafts_possible == []:
            return []

        #a = self.beads[index-2]
        b = self.beads[index-1]
        c = self.beads[index]
        #d = self.beads[index+1]
        e = self.beads[index+2]
        
        # try each of the four available pairs of crankshaft positions (including old one)
        # if both coords empty, then add them to filtered_crankshafts_possible
        filtered_crankshafts_possible = []
        for coord_pair in all_crankshafts_possible:
            if self.find_bead_at_this_coordinate(coord_pair[0]) == None and self.find_bead_at_this_coordinate(coord_pair[1]) == None:
                filtered_crankshafts_possible.append(coord_pair)

        # even though the position is possible, the *rotation* into that position may not be possible
        # b/c it's blocked, so this part checks for that
        final_possible_crankshafts = []
        for fil in filtered_crankshafts_possible:
            same_x, same_y, same_z, okay_to_copy = False, False, False, False
            if (fil[0][0] == c.x):
                same_x = True
            if (fil[0][1] == c.y):
                same_y = True
            if (fil[0][2] == c.z):
                same_z = True

            if (same_x and same_y) or (same_y and same_z) or (same_z and same_x):
                if (same_x):
                    okay_to_copy = okay_to_copy or (self.find_bead_at_this_coordinate((b.x+1, b.y, b.z)) == None and self.find_bead_at_this_coordinate((e.x+1, e.y, e.z)) == None)
                    okay_to_copy = okay_to_copy or (self.find_bead_at_this_coordinate((b.x-1, b.y, b.z)) == None and self.find_bead_at_this_coordinate((e.x-1, e.y, e.z)) == None)
                if (same_y):
                    okay_to_copy = okay_to_copy or (self.find_bead_at_this_coordinate((b.x, b.y+1, b.z)) == None and self.find_bead_at_this_coordinate((e.x, e.y+1, e.z)) == None)
                    okay_to_copy = okay_to_copy or (self.find_bead_at_this_coordinate((b.x, b.y-1, b.z)) == None and self.find_bead_at_this_coordinate((e.x, e.y-1, e.z)) == None)
                if (same_z):
                    okay_to_copy = okay_to_copy or (self.find_bead_at_this_coordinate((b.x, b.y, b.z+1)) == None and self.find_bead_at_this_coordinate((e.x, e.y, e.z+1)) == None)
                    okay_to_copy = okay_to_copy or (self.find_bead_at_this_coordinate((b.x, b.y, b.z-1)) == None and self.find_bead_at_this_coordinate((e.x, e.y, e.z-1)) == None)

            # if the crank position is NOT on opposite face, then turn is of course possible
            else:
                okay_to_copy = True

            # now we finally copy the pair to final_possible_crankshafts
            if (okay_to_copy):
                final_possible_crankshafts.append(fil)

        return final_possible_crankshafts

    def all_crankshafts_possible_for_bead(self, index):
        '''
            The crankshaft will be referred to as follows
                             (c = index)____________(d)
                                |                    |
                                |                    |
                                |                    |
                      (a)______(b)                  (e)______(f)

           This subproc finds out the common axis coords that (index-2), (index-1), (index+2), (index+3) have
           this subproc is to make crankshafts_possible look simpler
           THE ALGORITHM: if (index-1) and (index+2) are on the same Z coordinate (BUT (index) and (index+1)
           are not on that same Z coordinate), we try out combinations of 
           their neighbor coords for (index) and (index+1) as possible pairs for new crankshaft positions
           Repeat for X and Y axes
           ALTHOUGH COMPLEX, THIS CODE HAS BEEN TESTED         
        '''
        b = self.beads[index-1]
        c = self.beads[index]
        e = self.beads[index+2]
        
        code = 0
        if b.x  == e.x:
            code += 1
        if b.y  == e.y:
            code += 3
        if b.z  == e.z:
            code += 5

        # make sure the proposed "crank" is not linear, but an actual crank
        if code == 4 and b.x  == c.x and b.y  == c.y:
            return []
        if code == 6 and b.x  == c.x and b.z  == c.z:
            return []
        if code == 8 and b.y  == c.y and b.z  == c.z:
            return []

        possible_crankshaft_pairs = []
        if (code == 4 or code == 6):
            possible_crankshaft_pairs.append([(b.x+1, b.y, b.z), (e.x+1, e.y, e.z)])
            possible_crankshaft_pairs.append([(b.x-1, b.y, b.z), (e.x-1, e.y, e.z)])
        if (code == 4 or code == 8):
            possible_crankshaft_pairs.append([(b.x, b.y+1, b.z), (e.x, e.y+1, e.z)])
            possible_crankshaft_pairs.append([(b.x, b.y-1, b.z), (e.x, e.y-1, e.z)])
        if (code == 6 or code == 8):
            possible_crankshaft_pairs.append([(b.x, b.y, b.z+1), (e.x, e.y, e.z+1)])
            possible_crankshaft_pairs.append([(b.x, b.y, b.z-1), (e.x, e.y, e.z-1)])

        return possible_crankshaft_pairs

    def end_moves_possible_for_bead(self, index):      
        # in case function is called on a non-end protein, returns 0
        if index != 0 and index != self.N-1:
            return []
        
        end_moves_possible = []
        possible_coords = self.neighbor_coords_of_bead(self.beads[index])
        for pos in possible_coords:
            if self.find_bead_at_this_coordinate(pos) == None:
                end_moves_possible.append(pos)
        
        return end_moves_possible

    def choose_end_move(self, given_end_moves):
        random_index = int(self.rand.generate()*len(given_end_moves))
        return [given_end_moves[random_index]]
    
    def choose_crankshaft(self, given_crankshaft_moves):
        random_index = int(self.rand.generate()*len(given_crankshaft_moves))
        return given_crankshaft_moves[random_index]

    # returns a coordinate or pair of coordinates to move to
    def mc_move(self, index):
        possible_moves = self.end_moves_possible_for_bead(index)
        if possible_moves == []:
            possible_moves = self.corner_flip_possible_for_bead(index)
            if possible_moves == []:
                possible_move = self.crankshafts_possible_for_bead(index)
                if possible_moves == []:
                    return ([], '')
                else:
                    return (self.choose_crankshaft(possible_moves), "crankshaft")
            else:
                return (possible_moves, "corner flip")
        else:
            return (self.choose_end_move(possible_moves), "end move")

    def mc_accept(self, old_energy, new_energy):
        return self.rand.generate() < min(1.0, exp(-(new_energy - old_energy) / self.temperature))
