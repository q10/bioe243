//      HW#1 Question 2b.cpp
//      BE243
//      Stuart Howes <showes@berkeley.edu>

/*
rho* = rho/(m*sigma^3)

silver density = 1.537

*/

#include <iostream>
#include <valarray>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <cstdlib>
#include <math.h>
#include <cmath>
#include <string>

// using namespace std;


long MBIG = 1000000000;
long MSEED = time(NULL);
long MZ = 0;
float FAC = (1.0/MBIG);
float idum;

float RAN3() {
	// cout << idum << endl;
	static int inext,inextp;
	static float ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;
	if (idum < 0 || iff == 0) {
		iff=1;
		mj = MSEED-(idum < 0 ? -idum : idum);
		mj = mj % MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) {mk += MBIG;}
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) {ma[i] += MBIG;}
			}
		inext=0;
		inextp=31;
		idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) {mj += MBIG;}
	ma[inext]=mj;
	return mj*FAC;
	
}

int RANDINT (int low, int high) {
  return (int)(RAN3()*((float)high-low)) + low;
}

float min (float a, float b) {
  return (a < b ? a : b);
}

float BoxMueller() {
  static int turn = 0;
  static float x1, x2, r;
  
  if (turn == 1) {
    turn = 0;
    return x2*r;
  }
  else {
    r = 1.0;
    while (r >= 1.0) {
      x1 = (2.0*RAN3()) - 1.0;
      x2 = (2.0*RAN3()) - 1.0;
      r = x1*x1 + x2*x2;
    }
    r = sqrt(-2.0*log(r)/r);
    turn = 1;
    return x1*r;
  }
}


float force [14][108][3];				// Forces on each particle
float position [14][108][3];			// Position of current configuration of particles
float PositionAtTimeZero [108][3];			
float diffusitivity[650];
float diffusitivity_time[650];
float best_configuration [14][108][3];		// Position of best configurations of particles
float best_v [14];
float best_k [14];
float velocity [14][108][3];			// Velocities for all the particles in the current configuration
float rho_all [14][108];
float potential_energy [14];
float kinetic_energy [14];
float total_energy [14];
float current_temperature [14];
float setpoint_temperature [14];

float num_replicas = 1;
int N = 108;						// Number of particles for each configuration
int num_diffus_samples = 650;
float tau = 0.0001;                  // Parameter used with Bussi Thermostat, equilibration time


// PARTICLE PROPERTIES
float silver_e=2.5415E-3, silver_a=(4.09/3.20), silver_n=12, silver_m=6, silver_c=144.41, silver_density=1.53702761, silver_temperature=45.77;
float gold_e=1.2793E-2, gold_a=(4.08/2.70), gold_n=10, gold_m=8, gold_c=34.41, gold_density=0.92929022, gold_temperature=9.093;

void GetInitPositionsAndInitVelocities (int filecode) {
	float xi, yi, zi, Vx, Vy, Vz;
	int i = 0;

	FILE * filePosition;
	if (filecode==0)
	  filePosition = fopen("GOLD_BEST_CONFIG2.TXT", "r");
	else
	  filePosition = fopen("SILVER_BEST_CONFIG2.TXT", "r");
	if (filePosition == NULL) {
		std::cout << "Unable to open LJ_108_1.txt";
		exit(1); // terminate with error
        }
	fpos_t pos;
	while (!feof(filePosition)) {	
		fscanf (filePosition, "%f %f %f" , &xi, &yi, &zi);

                for (int h = 0; h < num_replicas; h++) {
			position[h][i][0] = xi;
			position[h][i][1] = yi;
			position[h][i][2] = zi;
			
			// save copy for best configuration
			best_configuration[h][i][0] = xi;
			best_configuration[h][i][1] = yi;
			best_configuration[h][i][2] = zi;

			float k = 1.0;
			//float k = sqrt (setpoint_temperature[h]);
			velocity[h][i][0] = k * BoxMueller();
			velocity[h][i][1] = k * BoxMueller();
			velocity[h][i][2] = k * BoxMueller();
		}
	
		i++;	
		fgetpos (filePosition, &pos);
		fscanf (filePosition, "%f %f %f" , &xi, &yi, &zi);
		if (!feof(filePosition)){ fsetpos (filePosition, &pos); }
		else break;
        }
	fclose (filePosition);
}


float GetRho(int particle_num, float box, float a, float m, int replica) {
  float rho = 0.0;

  for (int i = 0; i < N; i++){
    if (i != particle_num) {
      float deltaX = (position[replica][i][0] - position[replica][particle_num][0]);
      float deltaY = (position[replica][i][1] - position[replica][particle_num][1]);
      float deltaZ = (position[replica][i][2] - position[replica][particle_num][2]);

      //deltaX -= (box * round(deltaX/box));
      //deltaY -= (box * round(deltaY/box));
      //deltaZ -= (box * round(deltaZ/box));

      rho += pow(a / sqrt( pow(deltaX,2) + pow(deltaY,2) + pow(deltaZ,2) ), m);
    }
  }
  return rho;
}

void CalcForces(float box, float rcut, float a, float c, float m, float n){
  //Variables to calculate
  float r2 = 0, rcut2, Vcalc = 0, tmp_f, a_r = 0;
    
  // Clear forces from previous calculations
  // Use the same loop to calculate rho's for each and all particles
  for (int h=0; h<num_replicas; h++) {
    for (int i = 0; i < N; i++){
      for (int j = 0; j < 3; j++){
	force[h][i][j] = 0.0;
      }
      rho_all[h][i] = GetRho(i, box, a, m, h);
    }
  }

  for (int h=0; h<num_replicas; h++) {  
    potential_energy[h] = 0.0;

    for (int i = 0; i < N; i++){
      Vcalc = 0.0;

      for (int j = 0; j < N; j++){

	if (j != i){
	  float deltaX = (position[h][i][0] - position[h][j][0]);
	  float deltaY = (position[h][i][1] - position[h][j][1]);
	  float deltaZ = (position[h][i][2] - position[h][j][2]);
				
	  //deltaX -= (box * round(deltaX/box));
	  //deltaY -= (box * round(deltaY/box));
	  //deltaZ -= (box * round(deltaZ/box));
		
	  r2 = ( ( deltaX*deltaX ) +  ( deltaY*deltaY ) + ( deltaZ*deltaZ) );
	  //rcut2 = rcut * rcut;

	  //if (r2 < rcut2){
	  // Do calculation
	  a_r = a / sqrt(r2);
	  Vcalc += pow(a_r, n);
	    
	  tmp_f = ( n*pow(a_r, n) - (c/2.0)*m*(pow(rho_all[h][i],-0.5) + pow(rho_all[h][j],-0.5))*pow(a_r, m) ) / r2;

	  force[h][i][0] += deltaX * tmp_f;
	  force[h][i][1] += deltaY * tmp_f;
	  force[h][i][2] += deltaZ * tmp_f;
	  // remove double-counting
	  //force[h][j][0] -= deltaX * tmp_f;
	  //force[h][j][1] -= deltaY * tmp_f;
	  //force[h][j][2] -= deltaZ * tmp_f;
	  //}
	}
      }
      potential_energy[h] += 0.5*Vcalc - c * sqrt(rho_all[h][i]);
    }
  }

  /*
    // moved double-counting up
  for (int i = 0; i < N; i++){
    for (int j = 0; j < 3; j++){
      force[i][j] /= 2.0;
    }
  }
  */
}

void VelocityVerlet(float deltaT, float box, float rcut, float a, float c, float m, float n) {	
  float deltaT2 = deltaT * deltaT;
  
  for (int h=0; h<num_replicas; h++) {
    for (int i = 0; i < N; i++){
      // Update position
      position[h][i][0] += velocity[h][i][0]*deltaT + (1.0/2.0)*force[h][i][0]*deltaT2;
      position[h][i][1] += velocity[h][i][1]*deltaT + (1.0/2.0)*force[h][i][1]*deltaT2;
      position[h][i][2] += velocity[h][i][2]*deltaT + (1.0/2.0)*force[h][i][2]*deltaT2;
      
      // Update velocity at half step
      velocity[h][i][0] += force[h][i][0]*deltaT/2.0;
      velocity[h][i][1] += force[h][i][1]*deltaT/2.0;
      velocity[h][i][2] += force[h][i][2]*deltaT/2.0;

    }	
  }
  // New forces based on new positions
  CalcForces(box, rcut, a, c, m, n);	// Updates array of forces from F(t) to F(t + deltaT)
  
  // Update velocity at full step with new forces
  for (int h=0; h<num_replicas; h++) {
    for (int i = 0; i < N; i++){
      velocity[h][i][0] += force[h][i][0]*deltaT/2.0;
      velocity[h][i][1] += force[h][i][1]*deltaT/2.0;
      velocity[h][i][2] += force[h][i][2]*deltaT/2.0;
    }
  }
}

float GetKineticEnergy(int h) {
  float KineticEnergy = 0.0;
  for (int i = 0; i < N; i++){	 
    for (int j = 0; j < 3; j++){
      KineticEnergy += (1.0/2.0) * velocity[h][i][j] * velocity[h][i][j];
    }
  }
  return KineticEnergy;
}

void BussiThermostat(float deltaT) {
  /* Berendsen Thermostat
  // deltaT is the time step size
  float tau = 1.0;
  float lambda = sqrt(1 + deltaT/tau * (SetpointTemp/CurrentTemp - 1.0));
  for (int i = 0; i < N; i++){	 
  for (int j = 0; j < 3; j++){
  velocity[i][j] *= lambda;
  }
  }
  */
  
  float cc = exp(-deltaT/tau);

  for (int h=0; h<num_replicas; h++) {
    float d = (1-cc)*(setpoint_temperature[h]/current_temperature[h])/(N+1);
    float r = BoxMueller();
    float s = 0;
  
    // is it N-1, since it says from 1 to N-1 on the slides?
    for (int i=0; i<N-1; i++) {
      float si = BoxMueller();
      s += si*si;
    }
    float scale = sqrt( cc+(s+r*r)*d + 2.0*r*sqrt(cc*d) );
    //~ std::cout << "scale = " << scale << std::endl;

    if (r + sqrt(cc/d) < 0.0)
      scale = -scale;
  
    // rescale velocities
    for (int i = 0; i < N; i++){	 
      for (int j = 0; j < 3; j++){
	velocity[h][i][j] *= scale;
      }
    }
  }
}

float CalculateDiffusitivity(float timeEvolved, float box) {
  float diffusitivity = 0.0;
  
  for (int i=0; i<N; i++) {
    float deltaX = position[0][i][0] - PositionAtTimeZero[i][0];
    float deltaY = position[0][i][1] - PositionAtTimeZero[i][1];
    float deltaZ = position[0][i][2] - PositionAtTimeZero[i][2];

    //deltaX -= (box * round(deltaX/box));
    //deltaY -= (box * round(deltaY/box));
    //deltaZ -= (box * round(deltaZ/box));

    diffusitivity += deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;
  }

  diffusitivity /= (6*N*timeEvolved);
  return diffusitivity;
}

void update_energy_and_temperature () {
  for (int h=0; h < num_replicas; h++) {
    kinetic_energy[h] = GetKineticEnergy(h);									
    total_energy[h] = kinetic_energy[h] + potential_energy[h];
    current_temperature[h] = (2 * kinetic_energy[h])/((3 * N) - 6);
  }
}

void save_best_config () {
  for (int h=0; h < num_replicas; h++) {

    if (potential_energy[h] < best_v[h]){
      best_v[h] = potential_energy[h];
      best_k[h] = kinetic_energy[h];
      
      for (int v = 0; v < N; v++) {
	for (int w = 0; w < 3; w++) {
	  best_configuration[h][v][w] = position[h][v][w];
	}
      }
      
    }
  }
}

void setup_temperatures () {
  float temp_step = 0.3;
  float temperature = setpoint_temperature[0] - (temp_step*(num_replicas-2)/2.0);
  for (int h=1; h<num_replicas; h++) {
    setpoint_temperature[h] = temperature;
    temperature += temp_step;
  }
}

int get_best_global_config() {
  int best = 0;
  for (int h=0; h<num_replicas; h++) {
    if (best_v[h] < best_v[best])
      best = h;
  }
  return best;
}


void set_position_at_time_zero() {
  // save copy for diffusitivity calculation later
  for (int i=0; i<N; i++) {
    PositionAtTimeZero[i][0] = position[0][i][0];
    PositionAtTimeZero[i][1] = position[0][i][1];
    PositionAtTimeZero[i][2] = position[0][i][2];
  }
}

int main(int argc, char** argv){

  // declare metal properties
  float m, a, n, c, density;
  int filecode = 0;
  FILE * fileData;
  FILE * fileBestConfig;
  FILE * fileResults;
  FILE * fileDiffusion;
  
  // check for second and third arguments, and initialize metal properties accordingly
  if (argc < 2) {
    std::cout << "Program must be invoked with argument \'gold\' or \'silver\'\n\n";
    return -1; // terminate with error
  }
  else if (strcmp(argv[1],"gold") == 0) {
    m = gold_m;
    a = gold_a;
    n = gold_n;
    c = gold_c;
    density = gold_density;
    setpoint_temperature[0] = gold_temperature;
    //fileData = fopen("GOLD.RESULTS", "w");
    //fileBestConfig = fopen("GOLD_BEST_CONFIG.TXT", "w");
    //fileResults = fopen("GOLD.DATA", "w");
    fileDiffusion = fopen("GOLD_DIFFUSION.DATA2", "w");
  }
  else if (strcmp(argv[1],"silver") == 0) {
    m = silver_m;
    a = silver_a;
    n = silver_n;
    c = silver_c;
    density = silver_density;
    setpoint_temperature[0] = silver_temperature;
    //fileData = fopen("SILVER.RESULTS", "w");
    //fileBestConfig = fopen("SILVER_BEST_CONFIG.TXT", "w");
    //fileResults = fopen("SILVER.DATA", "w");
    fileDiffusion = fopen("SILVER_DIFFUSION.DATA2", "w");
    filecode=1;
  }
  else {
    std::cout << "Argument must be \'gold\' or \'silver\'\n\n";
    return -1;
  }

  // declaring variables, num timesteps to take, and d_s, the index of diffusitivity[]
  int timestart = time (NULL), timend, steps = 20000, d_s = 0;
  float timeEvolved, deltaT = 0.001;     // Size of time step
  float b = pow(abs(N/density) , (1.0/3.0));
  float rc = b/2;
  
  // check files are okay
  /*
  if (fileData == NULL) {
    std::cout << "Unable to open Results.txt" << std::endl;
    exit(1); // terminate with error
  }    
  if (fileBestConfig == NULL) {
    std::cout << "Unable to open BestConfig.txt" << std::endl;
    exit(1); // terminate with error
  } 
  if (fileResults == NULL) {
    std::cout << "Unable to open Data.txt" << std::endl;
    exit(1); // terminate with error
  }
  */
  if (fileDiffusion == NULL) {
    std::cout << "Unable to open Diffusion.txt" << std::endl;
    exit(1); // terminate with error
  }
    
  // Initialize positions, forces, energies, and temperatures
  setup_temperatures();
  GetInitPositionsAndInitVelocities(filecode);         // Start cold to prevent huge forces (should pass SetpointTemperature)
  CalcForces(b, rc, a, c, m, n);
  update_energy_and_temperature();
  
  // print results
  /*
  fprintf (fileResults, "NOTE: POTENTIAL ENERGY IS GIVEN BY THE SUTTON-CHEN POTENTIAL.  KINETIC ENERGY IS GIVEN BY THE PARTICLES' VELOCITIES.\n");
  fprintf (fileResults, "Box length = %f\n", b);
  fprintf (fileResults, "Starting Temperature = %f\n", current_temperature[0]);
  fprintf (fileData, "NOTE: POTENTIAL ENERGY IS GIVEN BY THE SUTTON-CHEN POTENTIAL.  KINETIC ENERGY IS GIVEN BY THE PARTICLES' VELOCITIES.\n");
  fprintf (fileData, "Time \t\tTotal Energy \tPotentialE(V) \tKineticE \tCurrent Temp \n");	// Print headers
  fprintf (fileData, "%f \t%f \t%f \t%f \t%f\n", 0.0, total_energy[0], potential_energy[0], kinetic_energy[0], current_temperature[0]);	// Print starting values

  */
  
  for (int i = 0; i < steps; i++){
    /*
    // keep the 0th config for MD, use the 1st to last configs as replicas
    if (RAN3() < 0.5) {
      int p = RANDINT(1, num_replicas);
      int q = p;
      if (p == 1)
	q++;
      else if (p == num_replicas-1 || RAN3() < 0.5)
	q--;
      else
	q++;

      float dBeta = (1.0/setpoint_temperature[q]) - (1.0/setpoint_temperature[p]);
      float dEnergy = potential_energy[q] - potential_energy[p];
      
      if (RAN3() < min(1.0, exp(dBeta*dEnergy))) {
	float temp;

	// swap configurations (positions, velocities, and forces)
	for (int y=0; y<N; y++) {
	  for (int z=0; z<3; z++) {
	    temp = position[p][y][z];
	    position[p][y][z] = position[q][y][z];
	    position[q][y][z] = temp;

	    temp = velocity[p][y][z];
	    velocity[p][y][z] = velocity[q][y][z];
	    velocity[q][y][z] = temp;

	    temp = force[p][y][z];
	    force[p][y][z] = force[q][y][z];
	    force[q][y][z] = temp;
	  }
	}
      }
    }
    */
    // Update position, force, velocity, energy, and temperatures
    VelocityVerlet(deltaT, b, rc, a, c, m, n);    
    update_energy_and_temperature();

    // Increase and control temperature, but wait until particles have spread equilibrated first
    if (i > 0.3*steps)
      BussiThermostat(deltaT);               

    // save best configuration
    save_best_config();
    
    timeEvolved = (i+1) * deltaT;      

    // set time zero for diffusitivity calculations, after equilibrium is reached (30000 steps)
    if (i == 0)
      set_position_at_time_zero();
    
    // calculate diffusitivity at intervals after equilibrium (30000 steps)
    if ((i%30 == 0) && (d_s < num_diffus_samples)) {
      diffusitivity[d_s] = CalculateDiffusitivity(timeEvolved, b);
      diffusitivity_time[d_s] = timeEvolved;
      d_s++;
    }

    //fprintf (fileData, "%f \t%f \t%f \t%f \t%f\n", timeEvolved, total_energy[0], potential_energy[0], kinetic_energy[0], current_temperature[0]);  
    
  }

  int BEST = get_best_global_config();
  /*  
  current_temperature[0] = (2 * kinetic_energy[0])/((3 * N) - 6);
  fprintf (fileResults, "Final Temperature = %f\n", current_temperature[0]);
  fprintf (fileResults, "Potential Energy = %f\n", potential_energy[0]);
  fprintf (fileResults, "Kinetic Energy = %f\n", kinetic_energy[0]);
  fprintf (fileResults, "Binding Energy (V/N) = %f\n", potential_energy[0] / N);
  //fprintf (fileResults, "Diffusitivity = %f\n\n", CalculateDiffusitivity(timeEvolved, b));

  fprintf (fileResults, "Global Minimum Search Results (found by Replica Exchange MD):\n");
  fprintf (fileResults, "Global Potential Energy Minimnum (V at best config) = %f\n", best_v[BEST]);
  fprintf (fileResults, "Kinetic Energy at Global Minimnum = %f\n", best_k[BEST]);
  fprintf (fileResults, "Total Energy at Global Minimnum = %f\n", best_v[BEST]+best_k[BEST]);
  fprintf (fileResults, "Binding Energy at Global Minimnum (V/N) = %f\n\n", best_v[BEST] / N);
    
  timend = time (NULL);
  fprintf (fileResults, "%i seconds to execute.", timend-timestart);
    
  fprintf (fileBestConfig, "Best Configuration:\n\n");
  fprintf (fileBestConfig, "x\t\ty\t\tz\n");
  for (int i = 0; i < N; i++) {
    fprintf (fileBestConfig, "%f \t%f \t%f\n", best_configuration[BEST][i][0], best_configuration[BEST][i][1], best_configuration[BEST][i][2]);
  }
  */
  // print diffusion data
  fprintf (fileDiffusion, "Diffusitivity Plot\n\n");
  fprintf (fileDiffusion, "Time \t\tDiffusitivity\n");
  for (int i = 0; i < num_diffus_samples; i++) {
    fprintf (fileDiffusion, "%f \t%f\n", diffusitivity_time[i], diffusitivity[i]);
  }

  //fclose (fileBestConfig);
  //fclose (fileData);
  //fclose (fileResults);
  fclose (fileDiffusion);
  
  return 0;
}

/*
  void calcDiff(int Switch, Coord* ct, int time, double dt, int N) {
    int t0index,i,j,n,CorrelTime;
    static int t0time[MAXT0],t0Counter,SampleCounter[MAXT];
    static double Vacf[MAXT],R2[MAXT];
    static Coord* c0 = new Coord[N*MAXT0];
    FILE *FilePtrMsd;
    switch(Switch) {
      // initialize everything
    case INITIALIZE:
      t0Counter=0;
      for(i=0;i<MAXT;i++) {
	  R2[i]=0.0;
	  SampleCounter[i]=0;
	}
      break;
    case SAMPLE:
      if((time % FREQT0)==0) {
	  // new time origin
	  // store the positions/velocities; the current velocities
	  // are ct.vx[i] and the current positions are ct.x[i] (and y and z).
	// question: why do you have to be careful with Pbc ?
	// 1. set t0index, t0counter and t0time[t0index]
	t0Counter++;
	t0index = (t0Counter - 1) % MAXT0;
	t0time[t0index] = time;
	// 2. store particle positions/velocities in c0, etc
	// to access position e.g particle i in x-dir, use ct[i].x
	// question: why should we use ct instead of c ?????
	for(i = 0; i < N; i++){
	  c0[i+N*t0index].x = ct[i].x;
	  c0[i+N*t0index].y = ct[i].y;
	  c0[i+N*t0index].z = ct[i].z;
	  // end modification
	}
    }
    // loop over all time origins that have been stored
    for(j=0;j<min(t0Counter,MAXT0);j++) {
      CorrelTime=time-t0time[j];
      // only if the time difference is shorter than the
      maximum correlation time
	// then add to R2
	if(CorrelTime < MAXT){
	  SampleCounter[CorrelTime]++;
	  for(n = 0; n < N; n++){
	    R2[CorrelTime] += (ct[n].x - c0[n+j*N].x)*(ct[n].x - c0[n+j*N].x);
	    R2[CorrelTime] += (ct[n].y - c0[n+j*N].y)*(ct[n].y - c0[n+j*N].y);
	    R2[CorrelTime] += (ct[n].z - c0[n+j*N].x)*(ct[n].z - c0[n+j*N].z);
	  }
	}
    }
    break;
  case WRITE_RESULTS:
    // write everything to disk
    FilePtrMsd=fopen("msd.dat","w");
    for(i=0;i<MAXT-1;i++) {
	if(SampleCounter[i]>0) {
	    R2[i]/=(double)(N*SampleCounter[i]);
	  }
	else{
	    R2[i]=0.0;
	}
	if(i>=(1.0/dt) )
	  fprintf(FilePtrMsd,"%lf %lf %lf %lf\n",(i+1)*dt,R2[i],R2[i]/(6.0*(i+1)*dt), (1.0/6.0) * (R2[i]-R2[(int)(1.0/dt)])/(dt*(i-(1.0/dt))) );
	else if(i>=1) {
	    fprintf(FilePtrMsd,"%lf %lf %lf\n",(i+1)*dt,R2[i],R2[i]/(6.0*(i+1)*dt));
	}
	else
	  fprintf(FilePtrMsd,"%lf %lf %lf\n",(i+1)*dt,R2[i],0.0);
    }
    fclose(FilePtrMsd);
    }
  }



 */
