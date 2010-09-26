#include "common.h"
#include "common.cpp"

#define MY_SEED 42709783

int main( int argc, char **argv ) {
  int nsteps=0, ntrials=100;
  
  if (argc < 2) {
    std::cout << "Program must be invoked with argument of form \'a\', \'b\', or \'c\'\n\n";
    return -1;
  }
  else if (strcmp(argv[1],"a") == 0) 
    nsteps=2;
  else if (strcmp(argv[1],"b") == 0) 
    nsteps=8;
  else if (strcmp(argv[1],"c") == 0)
    nsteps=32;
  else {
    std::cout << "Argument must be \'a\', \'b\', or \'c\'\n\n";
    return -1;
  }

  // set up and init all averages to zero
  double averages[ntrials];
  for (int i=0; i<ntrials; i++) {
    averages[i]=0.0;
  }

  // set up and init histogram to all zeros
  int histogram_intervals = 20;
  int histogram[histogram_intervals];
  for (int i=0; i<histogram_intervals; i++)
    histogram[i] = 0;

  // generate random numbers for nstep times
  for (int i=0; i<ntrials; i++) {
    for (int j=0; j<nsteps; j++) {
      averages[i] += ran3 (MY_SEED);
    }
  }
  
  // generate histogram of averages with 0.05-length intervals
  for (int i=0; i<ntrials; i++) {
    averages[i] /= (double)nsteps; // average the random numbers for each trial first
    for (int j=1; j<=histogram_intervals; j++) {
      if (averages[i] < (double)j/((double)histogram_intervals)) {
	histogram[j-1]++;
	break;
      }
    }
  }
  
  // print histogram
  for (int i=0; i<histogram_intervals; i++)
    std::cout << histogram[i] << std::endl;

  return 0;
}
