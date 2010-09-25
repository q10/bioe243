#include "common.h"
#include "common.cpp"

#define MY_SEED 42709783

std::vector<std::vector<int> > validPositions(int particle_positions[][2], int size, int current_x, int current_y);

int main( int argc, char **argv ) {
  bool type_one = true;

  if (argc < 2) {
    std::cout << "Program must be invoked with argument of form \'type1\' or \'type2\'\n\n";
    return -1;
  }
  else if (strcmp(argv[1],"type1") == 0) 
    type_one=true;
  else if (strcmp(argv[1],"type2") == 0) 
    type_one=false;
  else {
    std::cout << "Argument must be \'a\', \'b\', or \'c\'\n\n";
    return -1;
  }
  
  // if type 1 random walk, then we only need to save end position; else we need to save position on all steps
  int N = 100, particle_positions[N][2];
  bool failed_type_two = false;

  for (int conformation_number=0; conformation_number<100; conformation_number++) {
    
    // init array to zero to ensure clean slate (for easier debugging)
    for (int i=0; i<N; i++) 
      for (int j=0; j<2; j++)
	particle_positions[i][j] = 0;
    
    // grow protein, starting with array position 1 (position 0 is fixed to (0,0)
    for (int i=1; i<N; i++) {    
      // if we are running type-1 random walk
      if (type_one) {
	// set to old position first
	particle_positions[i][0] =  particle_positions[i-1][0];
	particle_positions[i][1] =  particle_positions[i-1][1];
	
	// equal probability of choosing X or Y axis, and of moving in + or - direction along that axis	
	int tmp_index = (ran3(MY_SEED)-0.5 < 0) ? 0 : 1;
	particle_positions[i][tmp_index] += (ran3(MY_SEED)-0.5 < 0) ? -1 : 1;
      }
      
      // else run type-2 random walk
      else {
	// get vector of valid next positions and choose one randomly
	std::vector<std::vector<int> > valid_positions = validPositions(particle_positions, i, particle_positions[i-1][0], particle_positions[i-1][1]);

	// if the protein curls unto itself and cannot random-walk anymore, then pull fail switch
	if (valid_positions.empty()) {
	  failed_type_two = true;
	  break;
	}
	
	// else choose next position
	else {
	  int rand_pos = (int)(ran3(MY_SEED)*valid_positions.size());
	  particle_positions[i][0] = valid_positions[rand_pos][0];
	  particle_positions[i][1] = valid_positions[rand_pos][1];
	}
      }
    }

    // if failed to produce full-length type-2 generated protein, then reset, and discount that conformation
    if (failed_type_two) {
      failed_type_two = false;
      conformation_number--;
      continue;
    }
    
    // after each conformation generated, print distant from end point to (0,0) (end-to-end distance)
    std::cout << setprecision(10)
	      << sqrt(pow((double)particle_positions[N-1][0],2) + 
		      pow((double)particle_positions[N-1][1],2)) << std::endl;
  }

  return 0;
}


std::vector<std::vector<int> > validPositions(int particle_positions[][2], int limited_size, int current_x, int current_y) {

  std::vector<std::vector<int> > valid_positions;
  
  int possible_positions[][2] = {{current_x-1, current_y},
				 {current_x+1, current_y},
				 {current_x, current_y-1},
				 {current_x, current_y+1}};
  bool possible[4];
  for (int i=0; i<4; i++) {
    possible[i] = true;
  }
  
  // check each of the 4 neighbors to see if available for next move
  for (int i=0; i<4; i++) {
    for (int j=0; j<limited_size; j++) {
      if ((possible_positions[i][0] == particle_positions[j][0]) && 
	  (possible_positions[i][1] == particle_positions[j][1])) {
	possible[i] = false;
	break;
      }
    }
  }
  
  // check the boolean flag to see which neighbor positions are available, and append those positions to vector
  for (int i=0; i<4; i++) {
    if (possible[i]) {
      std::vector<int> tmp_position;
      tmp_position.push_back(possible_positions[i][0]);
      tmp_position.push_back(possible_positions[i][1]);    
      valid_positions.push_back(tmp_position);
    }
  }
  return valid_positions;
}
