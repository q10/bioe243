//      HW#1 Question 2b.cpp
//      BE243
//      Stuart Howes <showes@berkeley.edu>

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

using namespace std;


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
	// cout << *idum;
	if (idum < 0 || iff == 0) {
		iff=1;
		mj = MSEED-(idum < 0 ? -idum : idum);
		//~ cout << "mj = " << mj << endl;
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

float density = 0.004511;		            // Density of fluid
float force [108][3];				// Forces on each particle
float position [108][3];			// Position of current configuration of particles
float PositionAtTimeZero [108][3];			
float AllPositions [1080][3];		// Position of all configurations of particles
float velocity [108][3];			// Velocities for all the particles in the current configuration
float velocityHalf [108][3];
float velocityNew [108][3];
float rho_all[108];
float SetpointTemp = 45; 
float Vfinal = 0;
int N = 108;						// Number of particles for each configuration
float tau = 0.001;                  // Parameter used with Bussi Thermostat, equilibration time

// PARTICLE PROPERTIES
float silver_e=2.5415E-3, silver_a=4.09, silver_n=12, silver_m=6, silver_c=144.41;
float gold_e=1.2793E-2, gold_a=4.08, gold_n=10, gold_m=8, gold_c=34.41;

int GetInitPositionsAndInitVelocities(){
	FILE * filePosition;
	//~ FILE * fileVelocity;
	float xi, yi, zi, Vx, Vy, Vz;
	int i = 0;
	
	filePosition = fopen("LJ_108_1.txt", "r");
	if (filePosition == NULL) {
		cout << "Unable to open LJ_108_1.txt";
		exit(1); // terminate with error
	}
	fpos_t pos;
	while (!feof(filePosition)) {	
			fscanf (filePosition, "%f %f %f" , &xi, &yi, &zi);
			position[i][0] = 5*xi;
			position[i][1] = 5*yi;
			position[i][2] = 5*zi;

			// save copy for diffusitivity calculation later
			PositionAtTimeZero[i][0] = xi;
			PositionAtTimeZero[i][1] = yi;
			PositionAtTimeZero[i][2] = zi;
			
			i++;	
			fgetpos (filePosition, &pos);
			fscanf (filePosition, "%f %f %f" , &xi, &yi, &zi);
			if (!feof(filePosition)){ fsetpos (filePosition, &pos); }
			else break;
	}
	fclose (filePosition);

	//~ i = 0;	
	//~ fileVelocity = fopen("init_velocities.out", "r");
	//~ if (fileVelocity == NULL) {
	        //~ cout << "Unable to open init_velocities.out"<<endl;
		//~ exit(1); // terminate with error
	//~ }
	//~ while (!feof(fileVelocity)) {	
			//~ fscanf (fileVelocity, "%f %f %f" , &Vx, &Vy, &Vz);
			//~ velocity[i][0] = Vx;
			//~ velocity[i][1] = Vy;
			//~ velocity[i][2] = Vz;
			//~ i++;	
			//~ fgetpos (fileVelocity, &pos);
			//~ fscanf (fileVelocity, "%f %f %f" , &Vx, &Vy, &Vz);
			//~ if (!feof(fileVelocity)) { fsetpos (fileVelocity, &pos); }
			//~ else break;
	//~ }
	//~ fclose (fileVelocity);	
    
    // Generate new velocities instead of reading in the old ones which were generated at a different temp.
    // This uses actual simulation temperature
    
    float u1, u2, u3, u4;			// Uniform random numbers to be pulled from RAN3
	float r;
	float g1, g2, g3, g4;			// Gaussian random number
	//~ float Vx, Vy, Vz;						// Velocities to be generated
   
    for (int i = 0; i < N; i++){		
        u1 = RAN3();
        u2 = RAN3();
        u3 = RAN3();
        u4 = RAN3();
        //~ u5 = RAN3();
        //~ u6 = RAN3();
        //~ printf ("Random values: u1 = %f, u2 = %f, u3 = %f, u4 = %f, u5 = %f and u6 = %f\n", u1, u2, u3, u4, u5, u6);
        
        u1 = (2.0 * u1) - 1.0;
        u2 = (2.0 * u2) - 1.0;
        u3 = (2.0 * u3) - 1.0;
        u4 = (2.0 * u4) - 1.0;
        //~ u5 = (2.0 * u5) - 1.0;
        //~ u6 = (2.0 * u6) - 1.0;
        //~ 
        //~ printf ("Random values transformed (-1, 1): u1 = %f, u2 = %f, u3 = %f, u4 = %f, u5 = %f and u6 = %f\n", u1, u2, u3, u4, u5, u6);
        
        first:
        r = (u1 * u1) + (u2 * u2);
        if ( (r < 1.0) && (r > 0.0) ){
            r = sqrt( -2.0*log(r)/r );
            g1 = u1 * r;
            g2 = u2 * r;
        }
        else {
            u1 = RAN3();
            u2 = RAN3();
            goto first;
        }
        
        second:
        r = (u3 * u3) + (u4 * u4);
        if ( (r < 1.0) && (r > 0.0) ){
            r = sqrt( -2.0*log(r)/r );
            g3 = u3 * r;
            g4 = u4 * r;
        }
        else {
            u3 = RAN3();
            u4 = RAN3();
            goto second;
        }
        
        //~ printf ("Gaussian values: g1 = %f, g2 = %f, g3 = %f, g4 = %f, g5 = %f and g6 = %f\n", g1, g2, g3, g4, g5, g6);
        
        // Debug...remove temp set
        SetpointTemp = 10;
        float k = sqrt (SetpointTemp);									// Generate actual velocities
            Vx = k * g1;
            Vy = k * g2;
            Vz = k * g3;
            
        //~ if ( u5 < 0 ) {Vx = -Vx;}									// Determine sign of velocities
        //~ if ( u6 < 0 ) {Vy = -Vy;}
        //~ if ( u7 < 0 ) {Vz = -Vz;}
        
        //~ printf ("Initial Velocities = %f, %f, %f \n", Vx, Vy, Vz);
        
        velocity[i][0] = Vx;
        velocity[i][1] = Vy;
        velocity[i][2] = Vz;
    }
}


float GetRho(int particle_num, float box, float a, float c, float m) {
  float rho = 0.0;

  for (int i = 0; i < N; i++){
    if (i != particle_num) {
      float deltaX = (position[i][0] - position[particle_num][0]);
      float deltaY = (position[i][1] - position[particle_num][1]);
      float deltaZ = (position[i][2] - position[particle_num][2]);

      deltaX -= (box * round(deltaX/box));
      deltaY -= (box * round(deltaY/box));
      deltaZ -= (box * round(deltaZ/box));

      rho += pow(a / sqrt( pow(deltaX,2) + pow(deltaY,2) + pow(deltaZ,2) ), m);
    }
  }
  return rho;
}

// THE NEW FORCE AND ENERGY FUNCTIONS
float CalcForces(float box, float rcut, float a, float c, float m, float n){
	//Variables to calculate
        float r2 = 0, r6 = 0, r12  = 0, rcut2, Vcalc = 0, Vcorrection = 0, tmp_f, a_r = 0;
    
	        // Clear forces from previous calculations
	        // Use the same loop to calculate rho's for each and all particles
		for (int i = 0; i < N; i++){
			for (int j = 0; j < 3; j++){
				force[i][j] = 0.0;
			}
			rho_all[i] = GetRho(i, box, a, c, m);
		}
	
		for (int i = 0; i < N; i++){
			for (int j = 0; j < N; j++){
            if (j != i){
				float deltaX = (position[i][0] - position[j][0]);
				float deltaY = (position[i][1] - position[j][1]);
				float deltaZ = (position[i][2] - position[j][2]);
				
				deltaX -= (box * round(deltaX/box));
				deltaY -= (box * round(deltaY/box));
				deltaZ -= (box * round(deltaZ/box));
		
				r2 = ( ( deltaX*deltaX ) +  ( deltaY*deltaY ) + ( deltaZ*deltaZ) );
				rcut2 = rcut * rcut;

				if (r2 < rcut2){
				  // Do calculation
				  a_r = a / sqrt(r2);
				  Vcalc += pow(a_r, n);
				  
				  // Reduced units, no e needed?
                  //tmp_f = e * ( n*pow(a_r, n) - c*m*(pow(GetRho(i),-0.5) + pow(GetRho(j),-0.5))*pow(a_r, m)/2.0 ) / r2;
				  tmp_f = ( n*pow(a_r, n) - (c/2.0)*m*(pow(rho_all[i],-0.5) + pow(rho_all[j],-0.5))*pow(a_r, m) ) / r2;

				  force[i][0] += deltaX * tmp_f;
				  force[i][1] += deltaY * tmp_f;
				  force[i][2] += deltaZ * tmp_f;
				  force[j][0] -= deltaX * tmp_f;
				  force[j][1] -= deltaY * tmp_f;
				  force[j][2] -= deltaZ * tmp_f;
				}
			}
            }
			Vcalc = 0.5*Vcalc - c * sqrt(rho_all[i]) ;
        }
		
		//~ float rcut3 = rcut * rcut * rcut;
		//~ float rcut9 = rcut3 * rcut3 * rcut3;
		//~ Vcorrection = 8 * 3.1415926535 * N * density * (1.0/(9 * rcut9) - 1.0/(3 * rcut3));
		//~ Vfinal = Vcalc + Vcorrection;
		Vfinal = Vcalc;
        
        
        
}

float VelocityVerlet(float deltaT, float box, float rcut, float a, float c, float m, float n){	
	float deltaT2 = deltaT * deltaT;
		
	for (int i = 0; i < N; i++){
		// Update position
		position[i][0] = position[i][0] + velocity[i][0]*deltaT + (1.0/2.0)*force[i][0]*deltaT2;
		position[i][1] = position[i][1] + velocity[i][1]*deltaT + (1.0/2.0)*force[i][1]*deltaT2;
		position[i][2] = position[i][2] + velocity[i][2]*deltaT + (1.0/2.0)*force[i][2]*deltaT2;

		// Update velocity at half step
		velocityHalf[i][0] = velocity[i][0] + force[i][0]*deltaT/2.0;
		velocityHalf[i][1] = velocity[i][1] + force[i][1]*deltaT/2.0;
		velocityHalf[i][2] = velocity[i][2] + force[i][2]*deltaT/2.0;

	}	
	
	// New forces based on new positions
	CalcForces(box, rcut, a, c, m, n);	// Updates array of forces from F(t) to F(t + deltaT)
		
	// Update velocity at full step with new forces
	for (int i = 0; i < N; i++){
		velocity[i][0] = velocityHalf[i][0] + force[i][0]*deltaT/2.0;
		velocity[i][1] = velocityHalf[i][1] + force[i][1]*deltaT/2.0;
		velocity[i][2] = velocityHalf[i][2] + force[i][2]*deltaT/2.0;
	}	
}

float GetKineticEnergy() {
  float KineticEnergy = 0.0;
  for (int i = 0; i < N; i++){	 
    for (int j = 0; j < 3; j++){
      KineticEnergy += (1.0/2.0) * velocity[i][j] * velocity[i][j];
    }
  }
  return KineticEnergy;
}

float BussiThermostat(float deltaT, float CurrentTemp) {
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
  float d = (1-cc)*(SetpointTemp/CurrentTemp)/(N+1);
  float r = BoxMueller();
  float s = 0;
  
  // is it N-1, since it says from 1 to N-1 on the slides?
  for (int i=0; i<N-1; i++) {
    float si = BoxMueller();
    s += si*si;
  }
  float scale = sqrt( cc+(s+r*r)*d + 2.0*r*sqrt(cc*d) );
  if (r + sqrt(cc/d) < 0.0)
    scale = -scale;

  // rescale velocities
  for (int i = 0; i < N; i++){	 
    for (int j = 0; j < 3; j++){
      velocity[i][j] *= scale;
    }
  }
}

float CalculateDiffusitivity(float deltaT, float step) {
  float diffusitivity = 0.0;
  
  for (int i=0; i<N; i++) {
    float deltaX = position[i][0] - PositionAtTimeZero[i][0];
    float deltaY = position[i][1] - PositionAtTimeZero[i][1];
    float deltaZ = position[i][2] - PositionAtTimeZero[i][2];
    diffusitivity += deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;
  }
  // average it and divide by 6*deltaT*step (also 6*total_time)
  diffusitivity /= (6*N*deltaT*step);
  return diffusitivity;
}

int main(int argc, char** argv){

  // declare metal properties
  float m, a, n, c;
  
  // check for second argument, and initialize metal properties accordingly
  if (argc < 2) {
    std::cout << "Program must be invoked with argument \'gold\' or \'silver\'" << std::endl;
    return -1; // terminate with error
  }
  else if (strcmp(argv[1],"gold") == 0) {
    m = gold_m;
    a = gold_a;
    n = gold_n;
    c = gold_c;
  }
  else if (strcmp(argv[1],"silver") == 0) {
    m = silver_m;
    a = silver_a;
    n = silver_n;
    c = silver_c;

  }
  else {
    std::cout << "Argument must be \'gold\' or \'silver\'\n\n";
    return -1;
  }

        // declaring variables
        int timestart = time (NULL), timend, steps = 10000;		                             // Number of time steps to take
		float KineticEnergy = 0, Temperature, TotalEnergy, timeEvolved, deltaT = 0.0001;     // Size of time step
		float b = pow(abs(N/density) , (1.0/3.0));
        float rc = b/2;
        
        printf("Box length = %f.\n", b);
        
		FILE * fileVelocityVerlet;
	 	fileVelocityVerlet = fopen("Results.txt", "w");
		if (fileVelocityVerlet == NULL) {
		  cout << "Unable to open VelocityVerlet.txt" << endl;
		  exit(1); // terminate with error
		}

		// Initialize positions and forces
		GetInitPositionsAndInitVelocities();
		CalcForces(b, rc, a, c, m, n);

        // tau = 
        
		// Calculate initial kinetic and total energy
		KineticEnergy = GetKineticEnergy();									
		TotalEnergy = KineticEnergy + Vfinal;
		Temperature = (2 * KineticEnergy)/(3 * N - 3);

		cout << "Starting Temperature = " << Temperature << endl;
		fprintf (fileVelocityVerlet, "Time, Total Energy, PotentialE(V), KineticE \n", TotalEnergy, Vfinal, KineticEnergy);	// Print starting values
        fprintf (fileVelocityVerlet, "0, %f, %f, %f\n", TotalEnergy, Vfinal, KineticEnergy);	// Print starting values
		
		for (int i = 0; i < steps; i++){
		        // Update position and velocity
			VelocityVerlet(deltaT, b, rc, a, c, m, n);						
                // Get Kinetic Energy
            KineticEnergy = GetKineticEnergy();
            Temperature = (2 * KineticEnergy)/(3 * N - 3);      // Calculate temperature
            BussiThermostat(deltaT, Temperature);               // Control temperature
			TotalEnergy = KineticEnergy + Vfinal;
			timeEvolved = i * deltaT;
            
            //~ printf ("t, TE, V, KE: %f, %f, %f, %f\n", timeEvolved, TotalEnergy, Vfinal, KineticEnergy);
			fprintf (fileVelocityVerlet, "%f, %f, %f, %f\n", timeEvolved, TotalEnergy, Vfinal, KineticEnergy);
		}
		
		Temperature = (2 * KineticEnergy)/(3 * N - 3);
		cout << "Final Temperature = " << Temperature << endl
		     << "Kinetic Energy = " << velocity[0][0] << endl
		     << "Binding Energy = " << TotalEnergy / N << endl
		     << "Diffusitivity = " << CalculateDiffusitivity(deltaT, steps) << endl;
		
	timend = time (NULL);
	cout << (timend-timestart) << " seconds to execute." << endl;
	
	//terminate the program
	return 0;
}
