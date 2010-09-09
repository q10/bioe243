#include "common.h"

#define BUF_SIZE 256
#define NUM_PARTICLES 108

int main( int argc, char **argv ) {

  // declare vars
  double particles[NUM_PARTICLES][3];
  double forces[NUM_PARTICLES][3];
  double LJEnergy;

  std::ifstream file;
  char line[BUF_SIZE];
  char *cToken;
  int i=0, j=0, NUM_CONFIGS;
  double tmp_x, tmp_y, tmp_z, dx, dy, dz, r, tmp_force;

  // attempt to open file
  if (argc < 2) {
    std::cout << "Program must be invoked with argument of form \'1\' or \'10\'\n" << std::endl;
    return -1;
  }
  
  if (atoi(argv[1]) == 1) {
    NUM_CONFIGS = 1;
    file.open( "LJ_108_1.txt" );
  }
  else if (atoi(argv[1]) == 10) {
    NUM_CONFIGS = 10;
    file.open( "LJ_108_10.txt" );
  }
  else {
    std::cout << "Argument must be \'1\' or \'10\'\n" << std::endl;
    return -1;
  }
    
  if( !file.good() ) {
    std::cout << "Failed to open \"LJ_108_1.txt / LJ_108_10.txt\" for input." << std::endl;
    return -1;
  }

  // run energy and force calculations
  for (int config_num=0; config_num<NUM_CONFIGS; config_num++) {

    // pipe numbers from text file into array 
    i = 0;
    while ( file.good() ) {
      if (i == NUM_PARTICLES) {
	break;
      }
      file.getline(line, BUF_SIZE);
      cToken = strtok( line, " \r\n" );
    
      if( cToken != NULL) {
	tmp_x = strtod(cToken, NULL);
	cToken = strtok( NULL, " \r\n" );
	tmp_y = strtod(cToken, NULL);
	cToken = strtok( NULL, " \r\n" );
	tmp_z = strtod(cToken, NULL);
      
	particles[i][0] = tmp_x;
	particles[i][1] = tmp_y;
	particles[i][2] = tmp_z;
	i++;
      }
    }

    // initialize energy and force vectors to zero
    LJEnergy = 0.0;
    for (i=0; i<NUM_PARTICLES; i++) {
      for (j=0; j<3; j++) {
	forces[i][j] = 0.0;
      }
    }  

    // calculate energy and forces
    for (i=0; i<NUM_PARTICLES-1; i++) {
      for (j=i+1; j<NUM_PARTICLES; j++) {

	dx = particles[i][0] - particles[j][0];
	dy = particles[i][1] - particles[j][1];
	dz = particles[i][2] - particles[j][2];
	r = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));

	LJEnergy += 4.0 * (pow((1/r),12) - pow((1/r),6));
	tmp_force = 24 * ((2 * pow((1/r),12)) - pow((1/r),6));
      
	forces[i][0] -= (tmp_force * dx / pow((1/r),2));
	forces[i][1] -= (tmp_force * dy / pow((1/r),2));
	forces[i][2] -= (tmp_force * dz / pow((1/r),2));

	forces[j][0] += (tmp_force * dx / pow((1/r),2));
	forces[j][1] += (tmp_force * dy / pow((1/r),2));
	forces[j][2] += (tmp_force * dz / pow((1/r),2));      
      }
    }

    // output numbers
    std::cout << "\nCONFIG #"
	      << config_num+1
	      << "\nLJ Energy = "
	      << setprecision(10) << LJEnergy
	      << "\nForce Vector on Particle 108 = ("
	      << particles[107][0] << ", " << particles[107][1] << ", " << particles[107][2] << ")\n\n";
  }  
  return 0;
}

