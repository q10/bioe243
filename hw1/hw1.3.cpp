#include "common.h"

#define BUF_SIZE 256
#define NUM_DRUG_MOLECS 1000000
#define MY_SEED 437197

double randu( unsigned int init );
double ran3( int init );

const unsigned int IMAX = (unsigned int)pow(2.0,31.0); 

int main( int argc, char **argv ) {

  // declare vars
  double RANDU_Particles[NUM_DRUG_MOLECS][3];
  double RAN3_Particles[NUM_DRUG_MOLECS][3];
  int circle_points=0, total_points=0;

  if (argc < 2) {
    std::cout << "Program must be invoked with argument of form \'randu\', \'ran3\', or \'circle\'\n\n";  
  }

  // generate random numbers with both RANDU and RAN3
  for (int i=0; i<NUM_DRUG_MOLECS; i++) {
    
    RANDU_Particles[i][0] = randu(MY_SEED);
    RANDU_Particles[i][1] = randu(MY_SEED);
    RANDU_Particles[i][2] = randu(MY_SEED);

    RAN3_Particles[i][0] = ran3(MY_SEED);
    RAN3_Particles[i][1] = ran3(MY_SEED);
    RAN3_Particles[i][2] = ran3(MY_SEED);

    // print values
    /*std::cout << setprecision(10)
	      << RANDU_Particles[i][0] << " "
	      << RANDU_Particles[i][1] << " " 
	      << RANDU_Particles[i][2] << "\t\t"
	      << RANDU_Particles[i][0] << " "
	      << RANDU_Particles[i][1] << " " 
	      << RANDU_Particles[i][2] << "\n";*/
  }

  double x, y, r;
  for (int i=0; i<NUM_DRUG_MOLECS; i++) {
    x = 2.0*ran3(MY_SEED) - 1.0;
    y = 2.0*ran3(MY_SEED) - 1.0;
    r = sqrt(pow(x,2) + pow(y,2));
    if (r < 1)
      circle_points++;
    total_points++;
  }
  std::cout << "Ratio of Circle to Total Points = " << (double)circle_points/(double)total_points << endl;

  return 0;
}

// The sequence of "random" numbers generated with seed = 1 agrees with existing data
double randu( unsigned int init ) {
  static bool started = false;
  static unsigned int iseed;
  
  // initialize random seed when called for the first time
  if (!started) {
    if (init % 2 == 0)
      init++;
    iseed = init;
    started = true;
  }

  iseed = (int)(((long long)65539 * (long long)iseed) % (long long)IMAX); // convert to long long for multiplying large numbers
  return (double)iseed / (double)IMAX;
}

// This code follows exactly the algorithm described in Numerical Recipes
double ran3( int init ) {
  static const int m = IMAX-1;
  static int mseed = 161803398;
  static int iseed = -1;
  static int rand_table[55];

  if (iseed < 0) {
    int mj = abs(mseed-abs(init)) % m;
    rand_table[54] = mj;
    int ii, mk = 1;

    // initialize rand_table using init as seed and large number mseed
    for (int i=0; i<54; i++) {
      int ii = (21*i) % 55;
      rand_table[ii] = mk;
      mk = mj - mk;
      if (mk < 0)
	mk += m;
      mj = rand_table[ii];
    }

    // randomize table warming up the generator
    for (int k=0; k<4; k++) {
      for (int i=0; i<55; i++) {
	rand_table[i] -= rand_table[1+((i+30)%55)];
	if (rand_table[i] < 0)
	  rand_table[i] += m;
      }
    }
    iseed = 0;
  }

  // generate the next random number
  rand_table[iseed] = (rand_table[iseed] - rand_table[(iseed+31) % 55]);
  if (rand_table[iseed] < 0)
    rand_table[iseed] += m;
  iseed = (iseed+1) % 55;

  return (double)rand_table[iseed] / (double)m;
}
