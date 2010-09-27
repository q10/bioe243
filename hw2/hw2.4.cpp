#include "common.h"
#include "common.cpp"

#define MY_SEED 279683937

// declare functions for use
void initialize_energies();
double energy_of_individual_amino_acid(int i);
void store_neighbors(int i);
int find_amino_acid_in_this_coordinate(int x, int y, int z);
void set_amino_acid_coords(int index, int new_x, int new_y, int new_z);
bool accepted(double new_U);
bool corner_flip_possible(int index);
int crankshafts_possible(int index);
void crankshafts_possible_subproc(int index);
int end_moves_possible(int index);
void flip_corner(int index);
void crankshaft(int index, int num_choices);
void make_end_move(int index, int num_choices);
void clear_old_data();

// global variables - very bad programming style, but code more optimized 
// (no need to pass-by-value or worry about pointer syntax)
const int NUM_AMINO_ACIDS = 36;
const char AMINO_ACID_TYPE[NUM_AMINO_ACIDS+1] = "HHPHHPPHPPPHPPHPHHHHPPPPHHPPHPHHHHPP";
int NCYCLES, amino_acids[NUM_AMINO_ACIDS][3], neighbors[6][3]; // 4th param of amino acids array is individual particle energy
double current_U, new_U, average_U, average_U_2, T, amino_acid_energies[NUM_AMINO_ACIDS];

// these variables used for storing next available coordinates for moving into, 
// and for storing the old coordinates in case move is rejected
// available_next_coords has 6 rows for 3 possible crankshafts
// current_coords_being_displaced has 2 possible amino acids moved max per MC move (crankshaft)
int current_coords_being_displaced[2][3];
int available_next_coords[6][3];
int crankshaft_possible_pairs[8][3];

int main( int argc, char **argv ) {

  
  
  initialize_energies();
  
  for (int i=0; i<NCYCLES; i++) {
    for (int j=0; j<NUM_AMINO_ACIDS; j++) {   
      
      double tmp_new_energy = (energy_of_individual_amino_acid(j) / 2.0);
      new_U = current_U - amino_acid_energies[j]  + tmp_new_energy;
      if (accepted(new_U)) {
	amino_acid_energies[j] = tmp_new_energy;
	current_U = new_U;
      }
      else {
	set_amino_acid_coords(j, current_coords_being_displaced[0][0], 
			      current_coords_being_displaced[0][1], 
			      current_coords_being_displaced[0][2]);
      }
      
      average_U += current_U;

      // for debugging purposes; can comment it out later to optimize runtime
      clear_old_data();
      
    }
  }
  return 0;
}



// initialize energies, coordinate vector, U (the energy of the system), and U_2 (U^2)
void initialize_energies() {
  current_U = 0.0, average_U = 0.0, average_U_2  = 0.0;
  for (int i=0; i<NUM_AMINO_ACIDS; i++) {    
    amino_acid_energies[i] = energy_of_individual_amino_acid(i);
    average_U += amino_acid_energies[i];
  }
  average_U /= 2.0; // account for double-counting the energy
  current_U = average_U;
  average_U_2 = pow(average_U, 2);

  return;
}

double energy_of_individual_amino_acid(int i) {
  double energy = 0.0;
  int index;
  // stores coordinates of immediate neighbors in array
  store_neighbors(i);
  
  // for each of the neigbor positions, find the index of amino acid that 
  // lives in it, if any, and add inter-amino-acid energy to total energy
  for (int j=0; j<6; j++) {
    index = find_amino_acid_in_this_coordinate(neighbors[j][0], neighbors[j][1], neighbors[j][2]);
    
    if (index >= 0) {
      if (AMINO_ACID_TYPE[i] == 'H')
	energy += (AMINO_ACID_TYPE[index] == 'H') ? -1.0 : 0.0;
      else
	energy += (AMINO_ACID_TYPE[index] == 'H') ? 0.0 : -0.5;
    }
    else
      energy += (AMINO_ACID_TYPE[i] == 'H') ? 0.5 : -0.5;
  }
  
  return energy;
}

// stores the coordinates of neigbor positions of amino acid with array index i 
void store_neighbors(int i) {
    neighbors[0][0] = amino_acids[i][0]+1;
    neighbors[0][1] = amino_acids[i][1];
    neighbors[0][2] = amino_acids[i][2];

    neighbors[1][0] = amino_acids[i][0]-1;
    neighbors[1][1] = amino_acids[i][1];
    neighbors[1][2] = amino_acids[i][2];

    neighbors[2][0] = amino_acids[i][0];
    neighbors[2][1] = amino_acids[i][1]+1;
    neighbors[2][2] = amino_acids[i][2];

    neighbors[3][0] = amino_acids[i][0];
    neighbors[3][1] = amino_acids[i][1]-1;
    neighbors[3][2] = amino_acids[i][2];

    neighbors[4][0] = amino_acids[i][0];
    neighbors[4][1] = amino_acids[i][1];
    neighbors[4][2] = amino_acids[i][2]+1;

    neighbors[5][0] = amino_acids[i][0];
    neighbors[5][1] = amino_acids[i][1];
    neighbors[5][2] = amino_acids[i][2]-1;

    return;
}

// finds the array index of the amino acid that lives in coordinates (x,y,z) or returns -1 otherwise
int find_amino_acid_in_this_coordinate(int x, int y, int z) {
  for (int i=0; i<NUM_AMINO_ACIDS; i++) {
    if (amino_acids[i][0] == x && amino_acids[i][1] == y && amino_acids[i][2] == z)
      return i;
  }
  return -1;
}

// searches the coordinates_containing_amino_acids vector for amino acid wth index index and gives it new coords
void set_amino_acid_coords(int index, int new_x, int new_y, int new_z) {
  amino_acids[index][0] = new_x;
  amino_acids[index][1] = new_y;
  amino_acids[index][2] = new_z;
  return;
}

// returns the acceptance of an MC move based on the potential energy of the system
bool accepted(double new_U) {
  return ran3(MY_SEED) < min(1.0, exp(-(new_U - current_U) / T));
}

// returns a true if a corner flip is possible, and 
// if true, places the other corner's coords into available_next_coords' row 0
bool corner_flip_possible(int index) {
  // find the opposite corner's coordinates
  if (amino_acids[index][0] == amino_acids[index+1][0] && 
      amino_acids[index][0] == amino_acids[index-1][0]) {
    if (amino_acids[index][1] == amino_acids[index-1][1]) {
      available_next_coords[0][1] = amino_acids[index+1][1];
      available_next_coords[0][2] = amino_acids[index-1][2];
    }
    else {
      available_next_coords[0][1] = amino_acids[index-1][1];
      available_next_coords[0][2] = amino_acids[index+1][2];
    }
  }

  else if (amino_acids[index][1] == amino_acids[index+1][1] && 
	   amino_acids[index][1] == amino_acids[index-1][1]) {    
    if (amino_acids[index][0] == amino_acids[index-1][0]) {
      available_next_coords[0][0] = amino_acids[index+1][0];
      available_next_coords[0][2] = amino_acids[index-1][2];
    }
    else {
      available_next_coords[0][0] = amino_acids[index-1][0];
      available_next_coords[0][2] = amino_acids[index+1][2];
    }
  }

  else {    
    if (amino_acids[index][0] == amino_acids[index-1][0]) {
      available_next_coords[0][0] = amino_acids[index+1][0];
      available_next_coords[0][1] = amino_acids[index-1][1];
    }
    else {
      available_next_coords[0][0] = amino_acids[index-1][0];
      available_next_coords[0][1] = amino_acids[index+1][1];
    }
  }

  // after we obtain the coordinates, find out if it's an empty space, or return false otherwise
  return find_amino_acid_in_this_coordinate(available_next_coords[0][0], 
					    available_next_coords[0][1], 
					    available_next_coords[0][2]) < 0;
}

// returns the number of crankshafts possible with amino acid index and index+1, and 
// places the new coords in pairs into available_next_coords' row 0-5
int crankshafts_possible(int index) {

  /* The crankshaft will be referred to as follows
                     (index)____________(index+1)
                        |                   |
                        |                   |
                        |                   |
     (index-2)______(index-1)           (index+2)______(index+3)
  */
  
  // The crankshaft cannot be performed on the first two or last two amino acids
  if (index < 2 || index > NUM_AMINO_ACIDS-3)
    return 0;
  
  int num_crankshafts_possible = 0;
  crankshafts_possible_subproc(index);
  for (int i=0; i<4; i++) {
    if (find_amino_acid_in_this_coordinate(crankshaft_possible_pairs[2*i][0], 
					   crankshaft_possible_pairs[2*i][1], 
					   crankshaft_possible_pairs[2*i][2]) < 0 &&
	find_amino_acid_in_this_coordinate(crankshaft_possible_pairs[(2*i)+1][0], 
					   crankshaft_possible_pairs[(2*i)+1][1], 
					   crankshaft_possible_pairs[(2*i)+1][2]) < 0) {
      // copy coords (index)
      available_next_coords[2*num_crankshafts_possible][0] = crankshaft_possible_pairs[2*i][0];
      available_next_coords[2*num_crankshafts_possible][1] = crankshaft_possible_pairs[2*i][1];
      available_next_coords[2*num_crankshafts_possible][2] = crankshaft_possible_pairs[2*i][2];

      // copy coords (index+1)
      available_next_coords[(2*num_crankshafts_possible)+1][0] = crankshaft_possible_pairs[(2*i)+1][0];
      available_next_coords[(2*num_crankshafts_possible)+1][1] = crankshaft_possible_pairs[(2*i)+1][1];
      available_next_coords[(2*num_crankshafts_possible)+1][2] = crankshaft_possible_pairs[(2*i)+1][2];
      
      // increment number of possible crankshafts
      num_crankshafts_possible++;
    }
  }
  return num_crankshafts_possible;
}

// this subproc finds out the common axis coords that (index-2), (index-1), (index+2), (index+3) have
// this subproc is to make crankshafts_possible look simpler
// if they have same X and Y, the code returned is 4
// if X and Z, code returned is 6
// if Y and Z, code returned is 8
void crankshafts_possible_subproc(int index) {
  int code = 0;
  if (amino_acids[index-2][0] == amino_acids[index-1][0] &&
      amino_acids[index-1][0]  == amino_acids[index+2][0] &&
      amino_acids[index+2][0] == amino_acids[index+3][0])
    code++;
  if (amino_acids[index-2][1] == amino_acids[index-1][1] &&
      amino_acids[index-1][1]  == amino_acids[index+2][1] &&
      amino_acids[index+2][1] == amino_acids[index+3][1])
    code += 3;
  if (amino_acids[index-2][2] == amino_acids[index-1][2] &&
      amino_acids[index-1][2]  == amino_acids[index+2][2] &&
      amino_acids[index+2][2] == amino_acids[index+3][2])
    code += 5;

  int i = 0;
  if (code == 4 || code == 6) {
    crankshaft_possible_pairs[i][0] = amino_acids[index-1][0]+1;
    crankshaft_possible_pairs[i][1] = amino_acids[index-1][1];
    crankshaft_possible_pairs[i++][2] = amino_acids[index-1][2];
    
    crankshaft_possible_pairs[i][0] = amino_acids[index+2][0]+1;
    crankshaft_possible_pairs[i][1] = amino_acids[index+2][1];
    crankshaft_possible_pairs[i++][2] = amino_acids[index+2][2];
  
    crankshaft_possible_pairs[i][0] = amino_acids[index-1][0]-1;
    crankshaft_possible_pairs[i][1] = amino_acids[index-1][1];
    crankshaft_possible_pairs[i++][2] = amino_acids[index-1][2];
    
    crankshaft_possible_pairs[i][0] = amino_acids[index+2][0]-1;
    crankshaft_possible_pairs[i][1] = amino_acids[index+2][1];
    crankshaft_possible_pairs[i++][2] = amino_acids[index+2][2];
  }
  if (code == 4 || code == 8) {
    crankshaft_possible_pairs[i][0] = amino_acids[index-1][0];
    crankshaft_possible_pairs[i][1] = amino_acids[index-1][1]+1;
    crankshaft_possible_pairs[i++][2] = amino_acids[index-1][2];
    
    crankshaft_possible_pairs[i][0] = amino_acids[index+2][0];
    crankshaft_possible_pairs[i][1] = amino_acids[index+2][1]+1;
    crankshaft_possible_pairs[i++][2] = amino_acids[index+2][2];
    
    crankshaft_possible_pairs[i][0] = amino_acids[index-1][0];
    crankshaft_possible_pairs[i][1] = amino_acids[index-1][1]-1;
    crankshaft_possible_pairs[i++][2] = amino_acids[index-1][2];
    
    crankshaft_possible_pairs[i][0] = amino_acids[index+2][0];
    crankshaft_possible_pairs[i][1] = amino_acids[index+2][1]-1;
    crankshaft_possible_pairs[i++][2] = amino_acids[index+2][2];
  }
  if (code == 6 || code == 8) {
    crankshaft_possible_pairs[i][0] = amino_acids[index-1][0];
    crankshaft_possible_pairs[i][1] = amino_acids[index-1][1];
    crankshaft_possible_pairs[i++][2] = amino_acids[index-1][2]+1;
    
    crankshaft_possible_pairs[i][0] = amino_acids[index+2][0];
    crankshaft_possible_pairs[i][1] = amino_acids[index+2][1];
    crankshaft_possible_pairs[i++][2] = amino_acids[index+2][2]+1;
    
    crankshaft_possible_pairs[i][0] = amino_acids[index-1][0];
    crankshaft_possible_pairs[i][1] = amino_acids[index-1][1];
    crankshaft_possible_pairs[i++][2] = amino_acids[index-1][2]-1;
    
    crankshaft_possible_pairs[i][0] = amino_acids[index+2][0];
    crankshaft_possible_pairs[i][1] = amino_acids[index+2][1];
    crankshaft_possible_pairs[i++][2] = amino_acids[index+2][2]-1;
  }
  return;
}

// returns the number of coords possible for end amino acid to move to and 
// places the new coords into available_next_coords
int end_moves_possible(int index) {
  int num_end_moves_possible = 0;

  // get neighbors
  store_neighbors(index);

  for (int j=0; j<6; j++) {
    // if there is no amino acid present in this coordinate, then store it as available
    if (find_amino_acid_in_this_coordinate(neighbors[j][0], neighbors[j][1], neighbors[j][2]) < 0) {
      available_next_coords[num_end_moves_possible][0] = neighbors[j][0];
      available_next_coords[num_end_moves_possible][1] = neighbors[j][1];
      available_next_coords[num_end_moves_possible][2] = neighbors[j][2];
      num_end_moves_possible++;
    }
  }

  return num_end_moves_possible;
}

// flip corner and save old coords into current_coords_being_displaced
void flip_corner(int index) {
  set_amino_acid_coords(index, available_next_coords[0][0], 
			available_next_coords[0][1], available_next_coords[0][2]);
  return;
}

// choose a crankshaft, perform crankshaft, and save old coords into current_coords_being_displaced
void crankshaft(int index, int num_choices) {
  int choice = (int)(ran3(MY_SEED) * num_choices);

  for (int k=0; k<2; k++) {
    current_coords_being_displaced[k][0] = amino_acids[index+k][0];
    current_coords_being_displaced[k][1] = amino_acids[index+k][1];
    current_coords_being_displaced[k][2] = amino_acids[index+k][2];

    set_amino_acid_coords(index+k, available_next_coords[(choice*2)+k][0], 
			  available_next_coords[(choice*2)+k][1], 
			  available_next_coords[(choice*2)+k][2]);
  }
  return;
}

// make end move and save old coords into current_coords_being_displaced
void make_end_move(int index, int num_choices) {
  int choice = (int)(ran3(MY_SEED) * num_choices);
  
  current_coords_being_displaced[0][0] = amino_acids[index][0];
  current_coords_being_displaced[0][1] = amino_acids[index][1];
  current_coords_being_displaced[0][2] = amino_acids[index][2];

  set_amino_acid_coords(index, available_next_coords[choice][0], 
			available_next_coords[choice][1], available_next_coords[choice][2]);
  return;
}

// clears current_coords_being_displaced, available_next_coords, neighbors, and new_U 
void clear_old_data() {
  for (int i=0; i<2; i++) {
    for (int j=0; j<3; j++)
      current_coords_being_displaced[i][j] = 0;
  }

  for (int i=0; i<6; i++) {
    for (int j=0; j<3; j++) {
      available_next_coords[i][j] = 0;
      neighbors[i][j] = 0;
    }
  }

  for (int i=0; i<8; i++)
    for (int j=0; j<8; j++)
    crankshaft_possible_pairs[i][j] = 0;

  new_U = 0.0;
  return;
}
