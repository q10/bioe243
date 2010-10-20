from math import sqrt, log, exp, pi
from common import Ran3, BoxMueller

RHO = 0.8442
DT = 0.001
DTT = DT**2.0
TARGET_TEMP = 0.742
NUM_PARTICLES = 108
# Because RHO* is 0.8442 and we assumed that sigma is 1,
# then RHO* = (N/V)(sigma^3) = N/V
# so V = N/RHO* = 127.93177
# so cube side length is 5.03878858

class Particle:
    def __init__(self, x, y, z, boundary):
        self.x, self.y, self.z = x, y, z
        self.vx, self.vy, self.vz = 0, 0, 0
        self.fx, self.fy, self.fz, = 0, 0, 0
        # mass, generalized to suit mass matrices
        self.mx, self.my, self.mz = 1, 1, 1
        self.energy = 0
        self.boundary = boundary
        self.enclose_to_boundary()

    def __str__(self):
        position = "Position (" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ")\n"
        velocity = "Velocity (" + str(self.vx) + ", " + str(self.vy) + ", " + str(self.vz) + ")\n"
        force = "Force ("+ str(self.fx) + ", " + str(self.fy) + ", " + str(self.fz) + ")\n"
        return position + velocity + force + "Energy = " + str(self.energy)

    def position(self):
        return (self.x, self.y, self.z)
    def velocity(self):
        return (self.vx, self.vy, self.vz)
    def force(self):
        return (self.fx, self.fy, self.fz)

    def enclose_to_boundary(self):
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
    
    # boundary conditions
    def add_position(self, dx, dy, dz):
        self.x += dx
        self.y += dy
        self.z += dz
        self.enclose_to_boundary()
    def add_velocity(self, dvx, dvy, dvz):
        self.vx += dvx
        self.vy += dvy
        self.vz += dvz
    def add_force(self, dfx, dfy, dfz):
        self.fx += dfx
        self.fy += dfy
        self.fz += dfz

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

        # initialize energies
        for first in self.particles:
            for other in self.particles:
                if other is not first:
                    self.add_energy(first, other)

        # initialize forces
        self.update_forces()
        
        # calculate energy truncation
        # the integral evaluates to (8/9)*PI*rho*N*((sigma/rc^9)-3(sigma/rc^3))
        r_cut = 1.0/self.cube_boundary
        self.energy_truncation = 8.0 * pi * len(self.particles) * RHO * ((r_cut**9.0)-(3.0*(r_cut**3.0))) / 9.0
        self.update_observables()

    def update_observables(self):
        self.update_system_potential_energy()
        self.update_system_kinetic_energy()
        self.system_total_energy = self.system_potential_energy + self.system_kinetic_energy + self.energy_truncation
        self.average_energy_per_particle = self.system_total_energy / len(self.particles)

    def velocity_rescale(self):
        # T* = Tk_b / e, so T* = sum(mv^2)/3N
        current_temp = 0.0
        for first in self.particles:
            current_temp += (first.mx*(first.vx**2)) + (first.my*(first.vy**2)) + (first.mz*(first.vz**2))
        current_temp /= 3*len(self.particles)
        
        if abs(current_temp-self.target_temp) > 0.001:
            gamma = sqrt(self.target_temp/current_temp)
            for first in self.particles:
                first.vx, first.vy, first.vz = first.vx*gamma, first.vy*gamma, first.vz*gamma
        
    def update_system_potential_energy(self):
        self.system_potential_energy = 0.0
        for first in self.particles:
            self.system_potential_energy += first.energy
        self.system_potential_energy /= 2.0

    def update_system_kinetic_energy(self):
        self.system_kinetic_energy = 0.0
        for first in self.particles:
            self.system_kinetic_energy += (first.mx*(first.vx**2.0)) + (first.my*(first.vy**2.0)) + (first.mz*(first.vz**2.0))
        self.system_kinetic_energy /= 2.0
        
    def add_energy(self, first, other):
        (dx, dy, dz) = self.mirror_convention(first.x-other.x, first.y-other.y, first.z-other.z)
        r = sqrt((dx**2.0)+(dy**2.0)+(dz**2.0))
        first.energy += 4.0 * (((1/r)**12.0) - ((1/r)**6.0))

    # calculate and add new force between two particles
    def add_force_between(self, first, other):        
        (dfx, dfy, dfz) = self.force_between(first, other)
        first.add_force(-dfx, -dfy, -dfz)
        other.add_force(dfx, dfy, dfz)

    # calculate force in between two particles        
    def force_between(self, first, other):
        (dx, dy, dz) = self.mirror_convention(first.x-other.x, first.y-other.y, first.z-other.z)
        r = sqrt((dx**2.0)+(dy**2.0)+(dz**2.0))
        tmp_force = -24.0 * ((2.0 * ((1/r)**14.0)) - ((1/r)**8.0))
        return (tmp_force*dx, tmp_force*dy, tmp_force*dz)

    # sets energies
    def update_energies(self):
        for first in self.particles:
            first.energy = 0.0
            for other in self.particles:
                if other is not first:
                    self.add_energy(first, other)

    # sets forces, assumes, they are pre-initialized to zero
    def update_forces(self):
        for tmp_j in range(len(self.particles)-1):
            for tmp_k in range(tmp_j+1, len(self.particles)):
                self.add_force_between(self.particles[tmp_j], self.particles[tmp_k])


    # Enforce mirror image convention
    def mirror_convention(self, dx, dy, dz):
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


class ForwardEuler(Algorithm):
    def __init__(self, filename="LJ_108_1.txt", temp=TARGET_TEMP, density=RHO, num=NUM_PARTICLES):
        Algorithm.__init__(self, filename, temp, density, num)

    def MD_step(self):
        self.time += DT
        for first in self.particles:            
            # update position
            first.add_position(first.vx*DT, first.vy*DT, first.vz*DT)

            # update velocity
            first.add_velocity(first.fx*DT/first.mx, first.fy*DT/first.my, first.fz*DT/first.mz)

        # update forces and energies
        self.update_forces()
        self.update_energies()

        #self.velocity_rescale()
        self.update_observables()


class VelocityVerlet(Algorithm):
    def __init__(self, filename="LJ_108_1.txt", temp=TARGET_TEMP, density=RHO, num=NUM_PARTICLES):
        Algorithm.__init__(self, filename, temp, density, num)
        
    def MD_step(self):
        self.time += DT
        for first in self.particles:           
            # second terms of Taylor expansion, for less computing in next 2 steps
            sx = first.fx/(2.0*first.mx)
            sy = first.fy/(2.0*first.my)
            sz = first.fz/(2.0*first.mz)
            
            # update position         
            dx = (first.vx*DT) + (sx*DTT)
            dy = (first.vy*DT) + (sy*DTT)
            dz = (first.vz*DT) + (sz*DTT)
            first.add_position(dx, dy, dz)

            # update velocity at half step with current forces
            dvx = sx*DT
            dvy = sy*DT
            dvz = sz*DT
            first.add_velocity(dvx, dvy, dvz)

            # reset force
            first.fx, first.fy, first.fz = 0.0, 0.0, 0.0
            
        # update forces and energies
        self.update_forces()
        self.update_energies()
        
        # update velocity at second half step with new forces
        for first in self.particles:
            dvx = (first.fx*DT)/(2*first.mx)
            dvy = (first.fy*DT)/(2*first.my)
            dvz = (first.fz*DT)/(2*first.mz)
            first.add_velocity(dvx, dvy, dvz)
            
        #self.velocity_rescale()
        self.update_observables()
