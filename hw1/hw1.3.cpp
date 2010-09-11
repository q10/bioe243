#include "common.h"

#define NUM_DRUG_MOLECS 1000000
#define MY_SEED 437197

double randu( unsigned int init );
double ran3( int init );

const unsigned int IMAX = (unsigned int)pow(2.0,31.0); 

int main( int argc, char **argv ) {
  if (argc < 2) {
    std::cout << "Program must be invoked with argument of form \'randu\', \'ran3\', or \'circle\'\n\n";
    return -1;
  }

  if ((strcmp(argv[1],"randu") == 0) || (strcmp(argv[1],"ran3") == 0)) {
    double tmp_x, tmp_y, tmp_z;

    // generate random numbers with either RANDU and RAN3
    for (int i=0; i<NUM_DRUG_MOLECS; i++) {

      if (strcmp(argv[1],"randu") == 0) {
	tmp_x = randu(MY_SEED);
	tmp_y = randu(MY_SEED);
	tmp_z = randu(MY_SEED);
      }
      else {
	tmp_x = ran3(MY_SEED);
	tmp_y = ran3(MY_SEED);
	tmp_z = ran3(MY_SEED);
      }
      
      // print values to stdout
      std::cout << setprecision(10)
		<< tmp_x << "\t"
		<< tmp_y << "\t" 
		<< tmp_z << "\n";
    }
  }

  else if (strcmp(argv[1],"circle") == 0) {   
    int circle_points=0, total_points=0;
    double x, y, r;

    for (int i=0; i<1000000000; i++) {
      x = 2.0*ran3(MY_SEED) - 1.0;
      y = 2.0*ran3(MY_SEED) - 1.0;
      r = sqrt(pow(x,2) + pow(y,2));

      if (r < 1)
	circle_points++;

      total_points++;
    
      if ((i==9) || (i==99) || (i==999) || (i==9999) ||
	  (i==99999) || (i==999999) || (i==9999999) || (i==99999999) ||
	  (i==999999999)) {
	std::cout << "Ratio of Circle to Total Points (PI/4) at "
		  << i+1 << " trials = "
		  << setprecision(10)
		  << (double)circle_points/(double)total_points << "\n";
      }
    }
  }

  else {
    std::cout << "Argument must be \'randu\', \'ran3\' or \'circle\'\n\n";
    return -1;
  }
  
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
