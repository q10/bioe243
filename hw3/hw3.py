from math import sqrt, log, exp, pi
from common import Ran3, BoxMueller

RHO = 0.8442
MY_SEED = 59375327573
DT = 0.001
# Because RHO* is 0.8442 and we assumed that sigma is 1,
# then RHO* = (N/V)(sigma^3) = N/V
# so V = N/RHO* = 127.93177
# so cube side length is 5.03878858
CUBE_BOUNDARY = 5.03878858
HALF_BOUNDARY = CUBE_BOUNDARY/2.0

class Particle:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z
        self.vx, self.vy, self.vz = 0, 0, 0
        self.fx, self.fy, self.fz, = 0, 0, 0        
        self.mx, self.my, self.mz = 1, 1, 1
        self.energy = 0

    def __str__(self):
        position = "Particle (" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ")\n"
        velocity = "Velocity (" + str(self.vx) + ", " + str(self.vy) + ", " + str(self.vz) + ")\n"
        force = "Force ("+ str(self.fx) + ", " + str(self.fy) + ", " + str(self.fz) + ")\n"
        return  position + velocity + force + "Energy = " + str(self.energy)
    
    # boundary conditions
    def add_position(self, dx, dy, dz):
        self.x += dx
        self.y += dy
        self.z += dz
        if self.x < -HALF_BOUNDARY:
            self.x %= HALF_BOUNDARY
        if self.x > HALF_BOUNDARY:
            self.x %= -HALF_BOUNDARY
        if self.y < -HALF_BOUNDARY:
            self.y %= HALF_BOUNDARY
        if self.y > HALF_BOUNDARY:
            self.y %= -HALF_BOUNDARY
        if self.z < -HALF_BOUNDARY:
            self.z %= HALF_BOUNDARY
        if self.z > HALF_BOUNDARY:
            self.z %= -HALF_BOUNDARY
    def add_velocity(self, dvx, dvy, dvz):
        self.vx += dvx
        self.vy += dvy
        self.vz += dvz
    def add_force(self, dfx, dfy, dfz):
        self.fx += dfx
        self.fy += dfy
        self.fz += dfz

class Algorithm:
    def __init__(self, filename):
        self.BoxMueller = BoxMueller(MY_SEED)
        self.system_energy, self.energy_truncation = 0, 0
        self.particles = []    

        # initialize particles and positions
        file = open(filename, 'r')
        line = file.readline()        
        while line != '':
            a = map(float, line.strip().replace(',', ' ').split())
            self.particles.append(Particle(a[0], a[1], a[2]))
            line = file.readline()

        # initialize forces and energies
        for first in self.particles:
            for other in self.particles:
                if other is not first:
                    self.add_force_between(first, other)
                    self.update_energy(first, other)

        # calculate energy truncation
        # the integral evaluates to (8/9)*PI*rho*N*((sigma/rc^9)-3(sigma/rc^3))
        r_cut = 1.0/HALF_BOUNDARY
        self.energy_truncation = 8.0 * pi * len(self.particles) * RHO * ((r_cut**9)-(3*(r_cut**3))) / 9.0
        
    def update_system_energy(self):
        self.system_energy = 0.0
        for first in self.particles:
            self.system_energy += first.energy
        self.system_energy /= 2.0

    def update_energy(self, first, other):
        (dx, dy, dz) = self.mirror_convention(first.x-other.x, first.y-other.y, first.z-other.z)
        r = sqrt((dx**2)+(dy**2)+(dz**2))        
        # calculate and add LJ potential
        first.energy += 4.0 * (((1/r)**12.0) - ((1/r)**6.0))

    def collect_old_df_values_for(self, first):
        old_df_values = []
        for other in self.particles:
            if other is not first:
                old_df_values.append(self.force_between(first, other))
        return old_df_values

    # calculate and add new force between two particles
    def add_force_between(self, first, other):        
        (dfx, dfy, dfz) = self.force_between(first, other)

        # account for double-counting of forces in the for-loop, by dividing force by 2
        dfx, dfy, dfz = dfx/2.0, dfy/2.0, dfz/2.0

        first.add_force(dfx, dfy, dfz)
        other.add_force(-dfx, -dfy, -dfz)

    # updates the force between two particles
    def update_force_between(self, first, other, old_df_values, i):
        # remove old force
        (dfx, dfy, dfz) = old_df_values[i]

        # account for double-counting of forces in the for-loop, by dividing force by 2
        dfx, dfy, dfz = dfx/2.0, dfy/2.0, dfz/2.0

        first.add_force(-dfx, -dfy, -dfz)
        other.add_force(dfx, dfy, dfz)

        # calculate and add new force
        self.add_force_between(first, other)

    # calculate force in between two particles        
    def force_between(self, first, other):
        (dx, dy, dz) = self.mirror_convention(first.x-other.x, first.y-other.y, first.z-other.z)
        r = sqrt((dx**2)+(dy**2)+(dz**2))
        tmp_force = 24 * ((2 * ((1/r)**14.0)) - ((1/r)**8.0))
        return (tmp_force*dx, tmp_force*dy, tmp_force*dz)

    # Enforce mirror image convention
    def mirror_convention(self, dx, dy, dz):
        if dx < -HALF_BOUNDARY:
            dx %= HALF_BOUNDARY
        if dx > HALF_BOUNDARY:
            dx %= -HALF_BOUNDARY
        if dy < -HALF_BOUNDARY:
            dy %= HALF_BOUNDARY
        if dy > HALF_BOUNDARY:
            dy %= -HALF_BOUNDARY
        if dz < -HALF_BOUNDARY:
            dz %= HALF_BOUNDARY
        if dz > HALF_BOUNDARY:
            dz %= -HALF_BOUNDARY
        return (dx, dy, dz)


class ForwardEuler(Algorithm):
    def __init__(self):
        Algorithm.__init__(self, filename)

    def MD_step(self):
        for first in self.particles:
            # collect old df values
            old_df_values = self.collect_old_df_values_for(first)
            
            # update position
            first.add_position(first.vx*DT, first.vy*DT, first.vz*DT)

            # update velocity
            first.add_velocity(first.fx*DT/first.mx, first.fy*DT/first.my, first.fz*DT/first.mz)            

            # erase old forces between this particle and the rest
            # replace with new forces
            # update energy of particle using same loop
            i, first.energy = 0, 0
            for other in self.particles:
                if other is not first:
                    self.update_force_between(first, other, old_df_values, i)
                    self.update_energy(first, other)
                    i+=1
   
    def averages(self):
        pass


class VelocityVerlet(Algorithm):
    def __init__(self, filename):
        Algorithm.__init__(self, filename)
        #initialize velocities
        
    def MD_step(self):
        DTT = DT**2
        for first in self.particles:
            # collect old df values
            old_df_values = self.collect_old_df_values_for(first)
            
            # second terms of Taylor expansion, for less computing in next 2 steps
            sx = first.fx/(2*first.mx)
            sy = first.fy/(2*first.my)
            sz = first.fz/(2*first.mz)
            
            # update position         
            dx = (first.vx*DT) + (sx*DTT)
            dy = (first.vy*DT) + (sy*DTT)
            dz = (first.vz*DT) + (sz*DTT)
            first.add_position(dx, dy, dz)

            # update velocity at half step with current forces
            dvx = sz*DT
            dvy = sz*DT
            dvz = sz*DT
            first.add_velocity(dvx, dvy, dvz)

            # erase old forces between this particle and the rest
            # use new positions to calculate new forces, and
            # use same loop to update energy of particle
            i, first.energy = 0, 0
            for other in self.particles:
                if other is not first:
                    self.update_force_between(first, other, old_df_values, i)
                    self.update_energy(first, other)
                    i+=1

            # update velocity at full step with new forces
            dvx = (first.fx*DT)/(2*first.mx)
            dvy = (first.fy*DT)/(2*first.my)
            dvz = (first.fz*DT)/(2*first.mz)
            first.add_velocity(dvx, dvy, dvz)

            # update total system energy
            self.update_system_energy()
    def averages(self):
        pass
