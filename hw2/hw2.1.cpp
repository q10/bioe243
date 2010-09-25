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

  double averages[ntrials];
  for (int i=0; i<ntrials; i++) {
    averages[i]=0.0;
  }

  for (int i=0; i<ntrials; i++) {
    for (int j=0; j<nsteps; j++) {
      averages[i] += ran3 (MY_SEED);
    }
  }

  for (int i=0; i<ntrials; i++) {
    std::cout << averages[i]/nsteps << std::endl;
  }

  return 0;
}
