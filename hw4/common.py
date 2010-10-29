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
        self.rho = density
        self.cube_boundary = ((self.N/self.rho)**(1.0/3.0))/2.0
        
        for particle in self.particles:
            particle.set_boundary(self.cube_boundary)

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
        #dx, dy, dz = first.x-other.x, first.y-other.y, first.z-other.z
        (dx, dy, dz) = self.__mirror_convention(first.x-other.x, first.y-other.y, first.z-other.z)
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
        #dx, dy, dz = first.x-other.x, first.y-other.y, first.z-other.z
        (dx, dy, dz) = self.__mirror_convention(first.x-other.x, first.y-other.y, first.z-other.z)
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




'''
class Algorithm:
    def __init__(self, filename, temp=TARGET_TEMP, density=RHO, num=NUM_PARTICLES):
        self.system_potential_energy, self.system_kinetic_energy, self.system_total_energy, self.energy_truncation = 0.0, 0.0, 0.0, 0.0
        self.time, self.system_average_energy_per_particle = 0.0, 0.0
        self.particles = []
        self.target_temp = temp
        self.N = num
        self.rho = density
        self.cube_boundary = ((self.N/self.rho)**(1.0/3.0))/2.0

        # initialize particles, positions, and velocities
        r_file = open(filename, 'r')
        v_file = open('init_velocities.out', 'r')
        r_line = r_file.readline()
        v_line = v_file.readline()

        j = 0
        while r_line != '' and v_line != '' and j < self.N:
            r = map(float, r_line.strip().replace(',', ' ').split())
            v = map(float, v_line.strip().replace(',', ' ').split())
            tmp = Particle(r[0], r[1], r[2], self.cube_boundary)
            tmp.add_velocity(v[0], v[1], v[2])
            self.particles.append(tmp)
            r_line = r_file.readline()
            v_line = v_file.readline()
            j+=1
        r_file.close()
        v_file.close()
'''
