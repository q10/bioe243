#include "common.h"

#define BUF_SIZE 256
#define NUM_AMINOACIDS 27

int main( int argc, char **argv ) {

  // declare vars
  int AminoAcids[NUM_AMINOACIDS][3];
  int energy;
  char AminoAcidType[NUM_AMINOACIDS+1] = "HPPHPHHPHPPPHPPHHHPPHHHHPPH";

  std::ifstream file;
  char line[BUF_SIZE];
  char *cToken;
  int i=0, j=0, tmp_x, tmp_y, tmp_z, dx, dy, dz, NUM_CONFIGS;

  // attempt to open file
  if (argc < 2) {
    std::cout << "Program must be invoked with argument of form \'1\' or \'10\'\n" << std::endl;
    return -1;
  }
  
  if (atoi(argv[1]) == 1) {
    NUM_CONFIGS = 1;
    file.open( "protein_lattice_config_1.txt" );
  }
  else if (atoi(argv[1]) == 10) {
    NUM_CONFIGS = 10;
    file.open( "protein_lattice_config_10.txt" );
  }
  else {
    std::cout << "Argument must be \'1\' or \'10\'\n" << std::endl;
    return -1;
  }
    
  if( !file.good() ) {
    std::cout << "Failed to open \"protein_lattice_config_1.txt / protein_lattice_config_10.txt\" for input." << std::endl;
    return -1;
  }

  // run energy and force calculations
  for (int config_num=0; config_num<NUM_CONFIGS; config_num++) {

    // pipe numbers from text file into array 
    i = 0;
    while ( file.good() ) {
      if (i == NUM_AMINOACIDS) {
	break;
      }
      file.getline(line, BUF_SIZE);
      cToken = strtok( line, " \r\n" );
      
      if( cToken != NULL) {
	tmp_x = atoi(cToken);
	cToken = strtok( NULL, " \r\n" );
	tmp_y = atoi(cToken);
	cToken = strtok( NULL, " \r\n" );
	tmp_z = atoi(cToken);
	
	AminoAcids[i][0] = tmp_x;
	AminoAcids[i][1] = tmp_y;
	AminoAcids[i][2] = tmp_z;
	i++;
      }
    }

    // init energy to zero
    energy = 0;
    
    // calculate energy
    for (i=0; i<NUM_AMINOACIDS-2; i++) {
      for (j=i+2; j<NUM_AMINOACIDS; j++) {

	dx = AminoAcids[i][0] - AminoAcids[j][0];
	dy = AminoAcids[i][1] - AminoAcids[j][1];
	dz = AminoAcids[i][2] - AminoAcids[j][2];

	// take absolute value without calling math library for speed
	if (dx < 0)
	  dx = -dx;
	if (dy < 0)
	  dy = -dy;
	if (dz < 0)
	  dz = -dz;

	if ((dx + dy + dz) < 2) {
	  if  ((AminoAcidType[i] == 'H' && AminoAcidType[j] == 'H') || (AminoAcidType[i] == 'P' && AminoAcidType[j] == 'P'))
	    energy -= 1;
	}
      }
    }

    // print value
    std::cout << "\nCONFIG #"
	      << config_num+1
	      << "\nEnergy = "
	      << setprecision(10) << energy << "\n";
  }
  return 0;
}

