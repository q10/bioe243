from math import sqrt, log, exp, pi
from common import Ran3, BoxMueller

DENSITY = 0.8442
DT = 0.001
DTT = DT**2.0
TARGET_TEMP = 0.742
NUM_PARTICLES = 108
SEED = 23823094327
SILVER_E, SILVER_A, SILVER_N, SILVER_M, SILVER_C, SILVER_DENSITY, SILVER_TEMP = 2.5415E-3, 4.09, 12, 6, 144.41, 0.153702761, 45.77
GOLD_E, GOLD_A, GOLD_N, GOLD_M, GOLD_C, GOLD_DENSITY, GOLD_TEMP = 1.2793E-2, 4.08, 10, 8, 34.41, 0.92929022, 9.093

# Because RHO* is 0.8442 and we assumed that sigma is 1,
# then RHO* = (N/V)(sigma^3) = N/V
# so V = N/RHO* = 127.93177
# so cube side length is 5.03878858

class Particle:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z
        self.vx, self.vy, self.vz = 0, 0, 0
        self.fx, self.fy, self.fz, = 0, 0, 0
        self.m = 1
        self.energy = 0
        self.x0, self.y0, self.z0 = self.x, self.y, self.z
        self.rho = 0

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
    
    def add_position(self, dx, dy, dz):
        self.x += dx
        self.y += dy
        self.z += dz
    def add_velocity(self, dvx, dvy, dvz):
        self.vx += dvx
        self.vy += dvy
        self.vz += dvz
    def add_force(self, dfx, dfy, dfz):
        self.fx += dfx
        self.fy += dfy
        self.fz += dfz

class SCCluster:
    def __init__(self, filename="LJ_108_1.txt", metal="gold"):
        self.system_potential_energy, self.system_kinetic_energy, self.system_total_energy, self.current_temperature = 0.0, 0.0, 0.0, 0.0
        self.time, self.system_average_energy_per_particle = 0.0, 0.0
        self.particles = []
        self.normal_dist = BoxMueller(SEED)
        self.N = 108
        self.tau = 0.001  # Bussi thermostat variable

        if metal == "silver":
            self.target_temp = SILVER_TEMP
            self.density = SILVER_DENSITY
            self.n = SILVER_N
            self.m = SILVER_M
            self.c = SILVER_C
        else:
            self.target_temp = GOLD_TEMP
            self.density = GOLD_DENSITY
            self.n = GOLD_N
            self.m = GOLD_M
            self.c = GOLD_C

        # actually half of the box's boundary
        self.cube_boundary = ((self.N/self.density)**(1.0/3.0))/2.0
        
        # initialize particles, positions, and velocities
        r_file = open(filename, 'r')
        r_line = r_file.readline()
        temp_sqrt = sqrt(self.target_temp)
        while r_line != '':
            r = map(float, r_line.strip().replace(',', ' ').split())
            tmp = Particle(r[0], r[1], r[2], self.cube_boundary)
            tmp.add_velocity(temp_sqrt*self.normal_dist.generate()[0], temp_sqrt*self.normal_dist.generate()[0], temp_sqrt*self.normal_dist.generate()[0])
            self.particles.append(tmp)
            r_line = r_file.readline()
        r_file.close()

        # initialize rho and forces
        self.update_forces_and_energies()
        
        self.update_observables()

        # for finding best global minimum energy and associated configuration
        self.best_config = []
        for p in self.particles:
                self.best_config.append(p.position())
        self.best_potential_energy = self.system_potential_energy
        self.best_total_energy = 0.0
        

    def update_observables(self):
        self.update_system_kinetic_energy()
        self.current_temperature = (2 * self.system_kinetic_energy) / (3*self.N - 3)
        self.system_total_energy = self.system_potential_energy + self.system_kinetic_energy
        self.average_energy_per_particle = self.system_total_energy / len(self.particles)

    def calculate_rho_for_particle(self, i):
        rho = 0.0;
        for j in range(self.N):
            if j != i:
                (dx, dy, dz) = self.mirror_convention(self.particles[j].x-self.particles[i].x, self.particles[j].y-self.particles[i].y, self.particles[j].z-self.particles[i].z)
                rho += (1 / sqrt((dx**2.0)+(dy**2.0)+(dz**2.0))) ** self.m;
        return rho;

    def Bussi_Thermostat(self):
        cc = exp(-DT/self.tau)
        d = (1-cc)*(self.target_temp/self.current_temperature)/(self.N+1)
        r = self.normal_dist.generate()[0]
        s = 0
  
        for i in range(self.N-1):
            si = self.normal_dist.generate()[0]
            s += si**2
          
        scale = sqrt( cc+(s+r*r)*d + 2.0*r*sqrt(cc*d) )
        if r + sqrt(cc/d) < 0.0:
            scale = -scale

        # rescale velocities
        for particle in self.particles:
            particle.vx *= scale
            particle.vy *= scale
            particle.vz *= scale

    def update_system_kinetic_energy(self):
        self.system_kinetic_energy = 0.0
        for first in self.particles:
            self.system_kinetic_energy += (first.m*(first.vx**2.0)) + (first.m*(first.vy**2.0)) + (first.m*(first.vz**2.0))
        self.system_kinetic_energy /= 2.0

    def update_forces_and_energies(self):
        self.system_potential_energy = 0.0
        
	# Clear forces from previous calculations
	# Use the same loop to calculate rho's for each and all particles
        for i in range(self.N):
            self.particles[i].fx = 0.0
            self.particles[i].fy = 0.0
            self.particles[i].fz = 0.0
            self.particles[i].rho = self.calculate_rho_for_particle(i)

        # calculate forces and energies	
        for i in range(self.N):
            for j in range(self.N):
                if j != i:
                    (dx, dy, dz) = self.mirror_convention(self.particles[j].x-self.particles[i].x, self.particles[j].y-self.particles[i].y, self.particles[j].z-self.particles[i].z)
                    r2 = (dx**2.0)+(dy**2.0)+(dz**2.0)
		
		    # Do calculation
                    r_i = 1 / sqrt(r2)
                    self.system_potential_energy += r_i**self.n
                    tmp_f = ( self.n*(r_i**self.n) - (self.c/2.0)*self.m*(self.particles[i].rho**-0.5 + self.particles[j].rho**-0.5)*(r_i**self.m) ) / r2

                    self.particles[i].fx += dx * tmp_f
                    self.particles[i].fy += dy * tmp_f
                    self.particles[i].fz += dz * tmp_f

            self.system_potential_energy /= 2
            self.system_potential_energy -= self.c * sqrt(self.particles[i].rho)

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

    def calculate_diffusitivity(self):
        diffusitivity = 0.0
        for i in range(self.N):
            dx = self.particles[i].x - self.particles[i].x0
            dy = self.particles[i].y - self.particles[i].y0
            dz = self.particles[i].z - self.particles[i].z0
            diffusitivity += dx**2 + dy**2 + dz**2
            #(dx, dy, dz) = self.mirror_convention(self.particles[j].x-self.particles[i].x, self.particles[j].y-self.particles[i].y, self.particles[j].z-self.particles[i].z)

            # average it and divide by 6*deltaT*step (also 6*total_time)
        diffusitivity /= (6*self.N*self.time)
        return diffusitivity

    def MD_step(self):
        self.time += DT
        for first in self.particles:           
            # second terms of Taylor expansion, for less computing in next 2 steps
            sx = first.fx/(2.0*first.m)
            sy = first.fy/(2.0*first.m)
            sz = first.fz/(2.0*first.m)
            
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
        self.update_forces_and_energies()

        # update velocity at second half step with new forces
        for first in self.particles:
            dvx = (first.fx*DT)/(2*first.m)
            dvy = (first.fy*DT)/(2*first.m)
            dvz = (first.fz*DT)/(2*first.m)
            first.add_velocity(dvx, dvy, dvz)

        self.update_observables()            
        self.Bussi_Thermostat()

        # save global minimum configuration
        if self.system_potential_energy < self.best_potential_energy:
            print "true"
            self.best_potential_energy = self.system_potential_energy
            self.best_total_energy = self.system_total_energy
            self.best_config = []
            for p in self.particles:
                self.best_config.append(p.position())
        
