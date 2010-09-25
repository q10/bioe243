#include "common.h"

const unsigned int IMAX = (unsigned int)pow(2.0,31.0);

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
      ii = (21*i) % 55;
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
