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

using namespace std;


long MBIG = 1000000000;
long MSEED = 574894209;
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

double density = 0.8442;			// Density of fluid
float force [108][3];				// Forces on each particle
float position [108][3];			// Position of current configuration of particles
//~ float positionNew [108][3];			
//~ float AllPositions [1080][3];		// Position of all configurations of particles
float velocity [108][3];			// Velocities for all the particles in the current configuration
float velocityHalf [108][3];
//~ float velocityNew [108][3];
float Temp = 0.742;
float Vfinal = 0;
float Kb = 1.3806503e-23;
int N = 108;						// Number of particles for each configuration
double a, rc;
// double box;

int GetInitPositionsAndInitVelocities(){
	FILE * filePosition;
	FILE * fileVelocity;
	float xi, yi, zi;
	float Vx, Vy, Vz;
	int i = 0;
	
	filePosition = fopen("/home/stuart/Documents/UC Berkeley/BE243/Homework/Homework 3/LJ_108_1.txt", "r");
	if (filePosition == NULL) {
		cout << "Unable to open LJ_108_1.txt";
		exit(1); // terminate with error
	}
	
	fpos_t pos;
	while (!feof(filePosition)) {	
			fscanf (filePosition, "%f %f %f" , &xi, &yi, &zi);
			//~ cout << xi << ", " << yi << ", " << zi << "." << endl;
			position[i][0] = xi;
			position[i][1] = yi;
			position[i][2] = zi;
			i++;	
			fgetpos (filePosition, &pos);
			fscanf (filePosition, "%f %f %f" , &xi, &yi, &zi);
			if (!feof(filePosition)){fsetpos (filePosition, &pos);}
			else break;
	}
	//~ cout << "Read in " << i << " position(s) successfully." << endl;
	
	//~ N = i;
	i = 0;
	fclose (filePosition);
	
	fileVelocity = fopen("/home/stuart/Documents/UC Berkeley/BE243/Homework/Homework 3/init_velocities.out", "r");
	if (filePosition == NULL) {
		cout << "Unable to open init_velocities.out";
		exit(1); // terminate with error
	}
	
	while (!feof(fileVelocity)) {	
			fscanf (fileVelocity, "%f %f %f" , &Vx, &Vy, &Vz);
			//~ cout << xi << ", " << yi << ", " << zi << "." << endl;
			velocity[i][0] = Vx;
			velocity[i][1] = Vy;
			velocity[i][2] = Vz;
			i++;	
			fgetpos (fileVelocity, &pos);
			fscanf (fileVelocity, "%f %f %f" , &Vx, &Vy, &Vz);
			if (!feof(fileVelocity)){fsetpos (fileVelocity, &pos);}
			else break;
	}
	//~ cout << "Read in " << i << " velocities(s) successfully." << endl;
	
	//~ N = i;
	i = 0;
	fclose (fileVelocity);
	
	//~ printf ("Initial: Particle %d has position %f, %f and %f\n",     10+1, position[10][0],     position[10][1],     position[10][2]);
	//~ printf ("Initial: Particle %d has velocity %f, %f and %f\n",     10+1, velocity[10][0],     velocity[10][1],     velocity[10][2]);
	//~ printf ("Initial: Particle %d has force %f, %f and %f\n",        10+1, force[10][0],        force[10][1],        force[10][2]);
	
	//~ for (int k = 0; k < N; k++){
		//~ printf ("Particle %d has position %f, %f and %f\n", k, position[k][0], position[k][1], position[k][2]);
		//~ printf ("Particle %d has velocity %f, %f and %f\n", k, velocity[k][0], velocity[k][1], velocity[k][2]);
	//~ }
	
}

float CalcForces(double box, double rcut){
	//Variables to calculate
		float r2 = 0;
		float r6 = 0;
		float r12  = 0;
		float rcut2;
		float Vcalc = 0;
		float Vcorrection = 0;
	
	// Clear forces from previous calculations
		for (int a = 0; a < N; a++){
			for (int b = 0; b < 3; b++){
				force[a][b] = 0.0;
			}
		}
	
		for (int i = 0; i < N-1; i++){
			for (int j = i+1; j < N; j++){
				float deltaX = (position[i][0] - position[j][0]);
				float deltaY = (position[i][1] - position[j][1]);
				float deltaZ = (position[i][2] - position[j][2]);
				
				//~ printf ("DeltaX = %f, DeltaY = %f, DeltaZ = %f.\n", deltaX, deltaY, deltaZ);
				deltaX = deltaX - (box * round(deltaX/box));
				deltaY = deltaY - (box * round(deltaY/box));
				deltaZ = deltaZ - (box * round(deltaZ/box));
				//~ printf ("Min Image Con:\n Box = %f\n DeltaX = %f, DeltaY = %f, DeltaZ = %f.\n", box, deltaX, deltaY, deltaZ);
		
				r2 = ( ( deltaX*deltaX ) +  ( deltaY*deltaY ) + ( deltaZ*deltaZ) );
				rcut2 = rcut * rcut;
				if (r2 < rcut2){
				//~ cout << "r2 = " << r2 << endl;
				r6 = r2 * r2 * r2;
				r12 = r6 * r6;
			
				// Do calculation
				//~ cout << " V current = " << V << endl;
				Vcalc = Vcalc + 4.0 * ( (1.0/r12) - (1.0/r6) ) ;
				
				
				//~ cout << "V updated = " << V << endl;
				force[i][0] = force[i][0] + ( 24.0*( (2.0*(1.0/r12)) - (1.0/r6) )*(deltaX/r2) );
				force[i][1] = force[i][1] + ( 24.0*( (2.0*(1.0/r12)) - (1.0/r6) )*(deltaY/r2) );
				force[i][2] = force[i][2] + ( 24.0*( (2.0*(1.0/r12)) - (1.0/r6) )*(deltaZ/r2) );
				force[j][0] = force[j][0] - ( 24.0*( (2.0*(1.0/r12)) - (1.0/r6) )*(deltaX/r2) );
				force[j][1] = force[j][1] - ( 24.0*( (2.0*(1.0/r12)) - (1.0/r6) )*(deltaY/r2) );
				force[j][2] = force[j][2] - ( 24.0*( (2.0*(1.0/r12)) - (1.0/r6) )*(deltaZ/r2) );
				// printf ("Forces on particle (%d) are X = %Lf, Y = %Lf and Z = %Lf.\n", (ci+1), ForceXi, ForceYi, ForceZi);
				}
			}
		
			//~ if (i == 10){ 		
				//~ printf ("Particle %d has position %f, %f and %f\n", i+1, position[i][0], position[i][1], position[i][2]);
				//~ printf ("Particle %d has velocity %f, %f and %f\n", i+1, velocity[i][0], velocity[i][1], velocity[i][2]);
			//~ }
		}
		
		float rcut3 = rcut * rcut * rcut;
		float rcut9 = rcut3 * rcut3 * rcut3;
		Vcorrection = 8 * 3.1415926535 * N * density * (1.0/(9 * rcut9) - 1.0/(3 * rcut3));
		Vfinal = Vcalc + Vcorrection;
		//~ cout << "Vcorrection = " << Vcorrection << endl;
		//~ printf("Potential = %f\n", Vfinal);
	
	
}

float ForwardEuler(float deltaT, double box, double rcut){
	
	//~ printf ("Start: Particle %d has position %f, %f and %f\n",     10+1, position[10][0],     position[10][1],     position[10][2]);
	//~ printf ("Start: Particle %d has velocity %f, %f and %f\n",     10+1, velocity[10][0],     velocity[10][1],     velocity[10][2]);
	//~ printf ("Start: Particle %d has force %f, %f and %f\n",        10+1, force[10][0],        force[10][1],        force[10][2]);
			
		// Update position
	for (int i = 0; i < N; i++){
		position[i][0] = position[i][0] + velocity[i][0]*deltaT;
		position[i][1] = position[i][1] + velocity[i][1]*deltaT;
		position[i][2] = position[i][2] + velocity[i][2]*deltaT;
	}
	
		// Update velocity
	for (int i = 0; i < N; i++){
		velocity[i][0] = velocity[i][0] + force[i][0]*deltaT;
		velocity[i][1] = velocity[i][1] + force[i][1]*deltaT;
		velocity[i][2] = velocity[i][2] + force[i][2]*deltaT;
	}
	
	//~ printf ("End: Particle %d has position %f, %f and %f\n",     10+1, position[10][0],     position[10][1],     position[10][2]);
	//~ printf ("End: Particle %d has velocity %f, %f and %f\n",     10+1, velocity[10][0],     velocity[10][1],     velocity[10][2]);
	//~ printf ("End: Particle %d has force %f, %f and %f\n",        10+1, force[10][0],        force[10][1],        force[10][2]);
}

float VelocityVerlet(float deltaT, double box, double rcut){	
	float deltaT2;
	deltaT2 = deltaT * deltaT;
	
	//~ printf ("Start: Particle %d has position %f, %f and %f\n",     10+1, position[10][0],     position[10][1],     position[10][2]);
	//~ printf ("Start: Particle %d has velocity %f, %f and %f\n",     10+1, velocity[10][0],     velocity[10][1],     velocity[10][2]);
	//~ printf ("Start: Particle %d has velocityHalf %f, %f and %f\n", 10+1, velocityHalf[10][0], velocityHalf[10][1], velocityHalf[10][2]);
	//~ printf ("Start: Particle %d has force %f, %f and %f\n",        10+1, force[10][0],        force[10][1],        force[10][2]);
	
		
	for (int i = 0; i < N; i++){
		// Update position
		position[i][0] = position[i][0] + velocity[i][0]*deltaT + (1.0/2.0)*force[i][0]*deltaT2;
		position[i][1] = position[i][1] + velocity[i][1]*deltaT + (1.0/2.0)*force[i][1]*deltaT2;
		position[i][2] = position[i][2] + velocity[i][2]*deltaT + (1.0/2.0)*force[i][2]*deltaT2;
	}
	
		// Update velocity at half step
	for (int i = 0; i < N; i++){
		velocityHalf[i][0] = velocity[i][0] + force[i][0]*deltaT/2.0;
		velocityHalf[i][1] = velocity[i][1] + force[i][1]*deltaT/2.0;
		velocityHalf[i][2] = velocity[i][2] + force[i][2]*deltaT/2.0;
	}
	
		// New forces based on new positions
	CalcForces(box, rcut);;	// Updates array of forces from F(t) to F(t + deltaT)
		
		// Update velocity at full step with new forces
	for (int i = 0; i < N; i++){
		velocity[i][0] = velocityHalf[i][0] + force[i][0]*deltaT/2.0;
		velocity[i][1] = velocityHalf[i][1] + force[i][1]*deltaT/2.0;
		velocity[i][2] = velocityHalf[i][2] + force[i][2]*deltaT/2.0;
	}	
	
	for (int k = 0; k < N; k++){
		printf ("Particle %d has position %f, %f and %f\n", k, position[k][0], position[k][1], position[k][2]);
		printf ("Particle %d has velocity %f, %f and %f\n", k, velocity[k][0], velocity[k][1], velocity[k][2]);
	}
	//~ printf ("End: Particle %d has position %f, %f and %f\n",     10+1, position[10][0],     position[10][1],     position[10][2]);
	//~ printf ("End: Particle %d has velocity %f, %f and %f\n",     10+1, velocity[10][0],     velocity[10][1],     velocity[10][2]);
	//~ printf ("End: Particle %d has velocityHalf %f, %f and %f\n", 10+1, velocityHalf[10][0], velocityHalf[10][1], velocityHalf[10][2]);
	//~ printf ("End: Particle %d has force %f, %f and %f\n",        10+1, force[10][0],        force[10][1],        force[10][2]);
}

int main(int argc, char** argv){
	// declaring variables
		int timestart = time (NULL);
		int timend;
		double b;                   // Box length for minimum image convention
		int steps = 50;				// Number of time steps to take
		float deltaT = 0.001;		    // Size of time step
		float KineticEnergy = 0;
		float Temperature;
		float TotalEnergy;
		
		a = (1.0/3.0);
		b = (pow (abs (N/density) , a));
		rc = b/2;
	 	
	 	FILE * fileForwardEuler;
	 	fileForwardEuler = fopen("/home/stuart/Documents/UC Berkeley/BE243/Homework/Homework 3/ForwardEuler.txt", "w");
		if (fileForwardEuler == NULL) {
			cout << "Unable to open ForwardEuler.txt";
			exit(1); // terminate with error
		}
	 		 	
		GetInitPositionsAndInitVelocities();
		CalcForces(b, rc);										// Populate forces for initial positions
		
		for (int n = 0; n < N; n++){							// Calculate initial kinetic energy
			for (int a = 0; a < 3; a++){
				KineticEnergy = KineticEnergy + (1.0/2.0) * velocity[n][a] * velocity[n][a];
			}
		}
		TotalEnergy = KineticEnergy + Vfinal;
		fprintf (fileForwardEuler, "0 %f %f %f\n", TotalEnergy, Vfinal, KineticEnergy);	// Print starting values
		Temperature = (2 * KineticEnergy)/(3 * N - 3);
		cout << "Starting Temperature = " << Temperature << endl;
		
		float timeEvolved;
		for (int i = 0; i < steps; i++){
			KineticEnergy = 0;									// Reset kinetic energy
			ForwardEuler(deltaT, b, rc);						// Update position and velocity
			CalcForces(b, rc);									// Calculate new forces
			for (int n = 0; n < N; n++){						// Calculate Kinetic Energy
				for (int a = 0; a < 3; a++){
					KineticEnergy = KineticEnergy + (1.0/2.0) * velocity[n][a] * velocity[n][a];
				}
			}
			timeEvolved = i * deltaT;
			TotalEnergy = KineticEnergy + Vfinal;
			fprintf (fileForwardEuler, "%f %f %f %f\n", timeEvolved, TotalEnergy, Vfinal, KineticEnergy);
			//~ cout << "Forward Euler Step " << i+1 << " complete." << endl;
		}
		Temperature = (2 * KineticEnergy)/(3 * N - 3);
		cout << "Final Temperature = " << Temperature << endl;
		
		// Reinitialize for next numerical integrator (big ol' reset!)
		FILE * fileVelocityVerlet;
	 	fileVelocityVerlet = fopen("/home/stuart/Documents/UC Berkeley/BE243/Homework/Homework 3/VelocityVerlet.txt", "w");
		if (fileVelocityVerlet == NULL) {
			cout << "Unable to open VelocityVerlet.txt";
			exit(1); // terminate with error
		}
		
		GetInitPositionsAndInitVelocities();	
		CalcForces(b, rc);										// Populate forces for initial positions
		KineticEnergy = 0;									
		for (int n = 0; n < N; n++){							// Calculate initial kinetic energy
			for (int a = 0; a < 3; a++){
				KineticEnergy = KineticEnergy + (1.0/2.0) * velocity[n][a] * velocity[n][a];
			}
		}
		TotalEnergy = KineticEnergy + Vfinal;
		fprintf (fileVelocityVerlet, "0 %f %f %f\n", TotalEnergy, Vfinal, KineticEnergy);	// Print starting values
		Temperature = (2 * KineticEnergy)/(3 * N - 3);
		cout << "Starting Temperature = " << Temperature << endl;
		
		for (int i = 0; i < steps; i++){
			KineticEnergy = 0;									// Reset kinetic energy
			VelocityVerlet(deltaT, b, rc);						// Update position and velocity
			for (int n = 0; n < N; n++){						// Calculate Kinetic Energy
				for (int a = 0; a < 3; a++){
					KineticEnergy = KineticEnergy + (1.0/2.0) * velocity[n][a] * velocity[n][a];
				}
			}
			timeEvolved = i * deltaT;
			TotalEnergy = KineticEnergy + Vfinal;
			fprintf (fileVelocityVerlet, "%f %f %f %f\n", timeEvolved, TotalEnergy, Vfinal, KineticEnergy);
			//~ cout << "Velocity Verlet Step " << i+1 << " complete." << endl;
		}
		
		Temperature = (2 * KineticEnergy)/(3 * N - 3);
		cout << "Final Temperature = " << Temperature << endl;
		
	timend = time (NULL);
	cout << (timend-timestart) << " seconds to execute (including waiting for user input).";
	
	//terminate the program
return 0;}
