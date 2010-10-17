#include "common.h"
#include "common.cpp"

#define MY_SEED 279683937
#define BUF_SIZE 256

// declare functions for use
void print_possible_moves();
int grab_input_data(int argc, char **argv);
void initialize_energies(int table_type);
bool mc_move(int index);
void mc_accept(int index, int table_type, int temperature, bool making_crankshaft_move);
void update_observables();
void calculate_and_print_values(int table_type, int temperature);
double energy_of_individual_amino_acid(int i, int table_type);
double energy_table_lookup_value(int i, int index, int table_type);
void store_neighbors(int i);
int find_amino_acid_in_this_coordinate(int x, int y, int z);
void set_amino_acid_coords(int index, int new_x, int new_y, int new_z);
bool accepted(double new_U, int temperature);
bool corner_flip_possible(int index);
int crankshafts_possible(int index);
int crankshafts_possible_subproc(int index);
int end_moves_possible(int index);
void flip_corner(int index);
void crankshaft(int index, int num_choices);
void make_end_move(int index, int num_choices);
void clear_old_data();

// global variables - very bad programming style, but code more optimized 
// (no need to pass-by-value or worry about pointer syntax)
const int NUM_AMINO_ACIDS = 36;
char AMINO_ACID_TYPE[NUM_AMINO_ACIDS+1] = "HHPHHPPHPPPHPPHPHHHHPPPPHHPPHPHHHHPP";
int NCYCLES=100000000, amino_acids[NUM_AMINO_ACIDS][3]; // 4th param of amino acids array is individual particle energy
int AMINO_ACIDS_TEMPLATE[NUM_AMINO_ACIDS][3]; // the original particle positions; for start of each simulation
double current_U, new_U, average_U, average_U_2, amino_acid_energies[NUM_AMINO_ACIDS];
double T[] = {0.6, 0.8, 0.85, 1.0, 1.2, 1.4, 1.6};

// these variables used for storing next available coordinates for moving into, 
// and for storing the old coordinates in case move is rejected
// available_next_coords has 5 rows for 5 possible end moves
// available_next_crankshaft_coords has 6 rows for 3 possible crankshafts
// current_coords_being_displaced has 2 possible amino acids moved max per MC move (crankshaft)
// neighbors store coordinates for the 6 neighbors
int current_coords_being_displaced[2][3], available_next_coords[5][3], available_next_crankshaft_coords[6][3], crankshaft_possible_pairs[8][3], neighbors[6][3];

int main( int argc, char **argv ) {

  // initialize system
  if (grab_input_data(argc, argv) < 0)
    return -1;
  
  for (int table_type=0; table_type<2; table_type++) {
    for (int temperature=0; temperature<7; temperature++) {

      initialize_energies(table_type);
      // print_possible_moves();
  
      // run simulation for NCYCLES
      bool making_crankshaft_move;

      for (int i=0; i<NCYCLES; i++) {
	for (int j=0; j<NUM_AMINO_ACIDS; j++) {   
	  // decide moves
	  making_crankshaft_move = mc_move(j);

	  // accept move
	  mc_accept(j, table_type, temperature, making_crankshaft_move);
      
	  // update temperature and such
	  update_observables();

	  // for debugging purposes; can comment it out later to optimize runtime
	  // clear_old_data();
	}
      }
  
      calculate_and_print_values(table_type, temperature);
    }
  }  
  
  return 0;
}

void calculate_and_print_values(int table_type, int temperature) {
  average_U /= (double)(NCYCLES * NUM_AMINO_ACIDS); 
  average_U_2 /= (double)(NCYCLES * NUM_AMINO_ACIDS);
  double C_T = (average_U_2 - pow(average_U, 2)) / (pow(T[temperature], 2));
  string str = (table_type < 1) ? "Solvent-less Model;  " : "Solvent Model;  ";
  cout << "Energy Model = " << str
       << "Temperature = " << T[temperature] 
       << ";  <U^2> = " << average_U_2 
       << ";  <U> = " << average_U
       << ";  C(T) = " << C_T << endl;
  return;
}

void print_possible_moves() {
  int count = 0;
  for (int j=0; j<NUM_AMINO_ACIDS; j++) {
    count = end_moves_possible(j);
    if (count > 0) {
      cout << "For bead #" << j+1 << " end move to ";
      for (int k=0; k<count; k++) {
	cout<< available_next_coords[k][0] << " "
	    << available_next_coords[k][1] << " "
	    << available_next_coords[k][2];
	if (k==count-1)
	  cout << endl;
	else
	  cout << ", ";
      }
    }

    bool test = corner_flip_possible(j);
    if (test) {
      cout << "For bead #" << j+1 << " corner flip to "
	   << available_next_coords[0][0] << " "
	   << available_next_coords[0][1] << " "
	   << available_next_coords[0][2] << endl;
      
    }
      
    count = crankshafts_possible(j);
    if (count > 0) {
      cout << "For bead #" << j+1 << " and " << j+2 << " : crankshaft move to: ";
      for (int k=0; k<count; k++) {
	cout << "("
	     << available_next_coords[2*k][0] << " "
	     << available_next_coords[2*k][1] << " "
	     << available_next_coords[2*k][2] << " and " 
	     << available_next_coords[(2*k)+1][0] << " "
	     << available_next_coords[(2*k)+1][1] << " "
	     << available_next_coords[(2*k)+1][2] << ")";
	if (k == count-1)
	  cout << endl;
	else
	  cout << " or ";
      }
    }
  }
  return;
}

int grab_input_data(int argc, char **argv) {
  std::ifstream file;
  
  // attempt to open file
  if (argc < 2) {
    std::cout << "Program must be invoked with argument of form \"input_file\"\n" << std::endl;
    return -1;
  }

  // open input file
  file.open(argv[1]);  
      
  if(!file.good()) {
    std::cout << "Failed to open \"" << argv[1] << "\" for input." << std::endl;
    return -1;
  }
  
  char line[BUF_SIZE],  *cToken;
  int tmp_x, tmp_y, tmp_z, i = 0;
  
  // To make the programming less complex, the file format will be as follows:
  // "3D coordinates"
  // The length of the protein chain, ncycles, and the bead type will be 
  // hard-coded into the program
  
  while ( file.good() ) {
    if (i == NUM_AMINO_ACIDS)
      break;
    file.getline(line, BUF_SIZE);  
    cToken = strtok( line, " \r\n" );

    if( cToken != NULL) {
      tmp_x = atoi(cToken);
      cToken = strtok( NULL, " \r\n" );
      tmp_y = atoi(cToken);
      cToken = strtok( NULL, " \r\n" );
      tmp_z = atoi(cToken);

      AMINO_ACIDS_TEMPLATE[i][0] = tmp_x;
      AMINO_ACIDS_TEMPLATE[i][1] = tmp_y;
      AMINO_ACIDS_TEMPLATE[i][2] = tmp_z;
      i++;
    }
  }
  
  return 0;
}


// initialize energies, coordinate vector, U (the energy of the system), and U_2 (U^2)
void initialize_energies(int table_type) {
  current_U = 0.0, average_U = 0.0, average_U_2  = 0.0;

  for (int i=0; i<NUM_AMINO_ACIDS; i++)
    for(int j=0; j<3; j++)
      amino_acids[i][j] = AMINO_ACIDS_TEMPLATE[i][j];

  for (int i=0; i<NUM_AMINO_ACIDS; i++) {    
    amino_acid_energies[i] = energy_of_individual_amino_acid(i, table_type);
    average_U += amino_acid_energies[i];
  }

  average_U /= 2.0; // account for double-counting the energy
  current_U = average_U;
  average_U_2 = pow(average_U, 2);

  return;
}

// makes mc move and returns true if a crankshaft move is made, b/c crankshafts are special cases
bool mc_move(int index) {
  bool making_crankshaft_move = false;
  
  // for end amino acid cases
  if (index == 0 || index == NUM_AMINO_ACIDS-1)
    make_end_move(index, end_moves_possible(index));

  // choose to flip corner
  else if (ran3(MY_SEED)*2 < 1.0) {
    if (corner_flip_possible(index))
      flip_corner(index);
  }
  
  // else attempt to make crankshaft
  else {
    crankshaft(index, crankshafts_possible(index));
    making_crankshaft_move = true;
  }
  
  return making_crankshaft_move;
}

// accept move
void mc_accept(int index, int table_type, int temperature, bool making_crankshaft_move) {

  // calculate energy difference
  double tmp_new_energy_2, tmp_new_energy = (energy_of_individual_amino_acid(index, table_type) / 2.0);
  new_U = current_U - amino_acid_energies[index] + tmp_new_energy;
  if (making_crankshaft_move) {
    tmp_new_energy_2 = (energy_of_individual_amino_acid(index+1, table_type) / 2.0);
    new_U = new_U - amino_acid_energies[index+1] + tmp_new_energy_2;
  }
  
  // if move is accepted, then update the individual and total energies
  if (accepted(new_U, temperature)) {
    amino_acid_energies[index] = tmp_new_energy;
    if (making_crankshaft_move)
      amino_acid_energies[index+1] = tmp_new_energy_2;
    current_U = new_U;
  }

  // otherwise, reset amino coords to the old ones
  else {
    set_amino_acid_coords(index, current_coords_being_displaced[0][0], 
			  current_coords_being_displaced[0][1], 
			  current_coords_being_displaced[0][2]);
    if (making_crankshaft_move)
      set_amino_acid_coords(index+1, current_coords_being_displaced[1][0], 
			    current_coords_being_displaced[1][1], 
			    current_coords_being_displaced[1][2]);
	
  }
  return;
}

void update_observables() {
  average_U += current_U;
  average_U_2 += pow(current_U, 2);
  return;
}

double energy_of_individual_amino_acid(int i, int table_type) {
  double energy = 0.0;
  int index;
  // stores coordinates of immediate neighbors in array
  store_neighbors(i);
  
  // for each of the neigbor positions, find the index of amino acid that 
  // lives in it, if any, and add inter-amino-acid energy to total energy
  for (int j=0; j<6; j++) {
    index = find_amino_acid_in_this_coordinate(neighbors[j][0], neighbors[j][1], neighbors[j][2]);
    energy += energy_table_lookup_value(j, index, table_type);
  }
  
  return energy;
}

double energy_table_lookup_value(int i, int index, int table_type) {
  double energy = 0.0;

  // if using energy model 0 (w/o solvent) and index is legit
  if (table_type == 0 && index >= 0) {
    if (AMINO_ACID_TYPE[i] == 'H')
      energy += (AMINO_ACID_TYPE[index] == 'H') ? -1.0 : 0.0;
    else
      energy += (AMINO_ACID_TYPE[index] == 'P') ? -1.0 : 0.0;
  }
  
  // else we are using energy model 1
  else {
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
bool accepted(double new_U, int temperature) {
  return ran3(MY_SEED) < min(1.0, exp(-(new_U - current_U) / T[temperature]));
}

// returns a true if a corner flip is possible, and 
// if true, places the other corner's coords into available_next_coords' row 0
bool corner_flip_possible(int index) {
  
  // account for end molecules
  if (index == 0 || index == NUM_AMINO_ACIDS-1)
    return 0;
  
  // find the opposite corner's coordinates
  if (amino_acids[index][0] == amino_acids[index+1][0] && 
      amino_acids[index][0] == amino_acids[index-1][0]) {
    available_next_coords[0][0] = amino_acids[index][0];
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
    available_next_coords[0][1] = amino_acids[index][1];
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
    available_next_coords[0][2] = amino_acids[index][2];
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

// returns the number of crankshafts possible with amino acid index and index+1 
// (just the ones that are both open and possible to rotate to), and 
// places the new coords in pairs into available_next_coords' row 0-5
int crankshafts_possible(int index) {

  /* The crankshaft will be referred to as follows
                     (index)____________(index+1)
                        |                   |
                        |                   |
                        |                   |
     (index-2)______(index-1)           (index+2)______(index+3)
  */
  
  // The crankshaft cannot be performed on the first or last two amino acids
  if (index < 2 || index > NUM_AMINO_ACIDS-3)
    return 0;
  
  int num_trial_crankshafts_possible = crankshafts_possible_subproc(index), num_valid_crankshafts_possible = 0;
  bool same_x, same_y, same_z, okay_to_copy;

  // try each of the four available pairs of crankshaft positions (including old one)
  // if both coords empty, then add them to available_next_coords
  for (int i=0; i<num_trial_crankshafts_possible; i++) {
    if (find_amino_acid_in_this_coordinate(crankshaft_possible_pairs[2*i][0], 
					   crankshaft_possible_pairs[2*i][1], 
					   crankshaft_possible_pairs[2*i][2]) < 0 &&
	find_amino_acid_in_this_coordinate(crankshaft_possible_pairs[(2*i)+1][0], 
					   crankshaft_possible_pairs[(2*i)+1][1], 
					   crankshaft_possible_pairs[(2*i)+1][2]) < 0) {

      // even though the position is possible, the *rotation* into that position may not be
      // so this monstrous if statement checks for that
      same_x = false, same_y = false, same_z = false, okay_to_copy = false;
      if (crankshaft_possible_pairs[2*i][0] == amino_acids[index][0])
	same_x = true;
      if (crankshaft_possible_pairs[2*i][1] == amino_acids[index][1])
	same_y = true;
      if (crankshaft_possible_pairs[2*i][2] == amino_acids[index][2])
	same_z = true;

      if ((same_x && same_y) || (same_y && same_z) || (same_z && same_x)) {
	if (same_x)
	  okay_to_copy = okay_to_copy || (find_amino_acid_in_this_coordinate(amino_acids[index-1][0]+1, 
									     amino_acids[index-1][1], 
									     amino_acids[index-1][2]) < 0 &&
					  find_amino_acid_in_this_coordinate(amino_acids[index+2][0]+1, 
									     amino_acids[index+2][1], 
									     amino_acids[index+2][2]) < 0) || 
	    (find_amino_acid_in_this_coordinate(amino_acids[index-1][0]-1, 
						amino_acids[index-1][1], 
						amino_acids[index-1][2]) < 0 &&
	     find_amino_acid_in_this_coordinate(amino_acids[index+2][0]-1, 
						amino_acids[index+2][1], 
						amino_acids[index+2][2]) < 0);
	if (same_y)
	  okay_to_copy = okay_to_copy || (find_amino_acid_in_this_coordinate(amino_acids[index-1][0], 
									     amino_acids[index-1][1]+1, 
									     amino_acids[index-1][2]) < 0 &&
					  find_amino_acid_in_this_coordinate(amino_acids[index+2][0], 
									     amino_acids[index+2][1]+1, 
									     amino_acids[index+2][2]) < 0) || 
	    (find_amino_acid_in_this_coordinate(amino_acids[index-1][0], 
						amino_acids[index-1][1]-1, 
						amino_acids[index-1][2]) < 0 &&
	     find_amino_acid_in_this_coordinate(amino_acids[index+2][0], 
						amino_acids[index+2][1]-1, 
						amino_acids[index+2][2]) < 0);
	if (same_z)
	  okay_to_copy = okay_to_copy || (find_amino_acid_in_this_coordinate(amino_acids[index-1][0], 
									     amino_acids[index-1][1], 
									     amino_acids[index-1][2]+1) < 0 &&
					  find_amino_acid_in_this_coordinate(amino_acids[index+2][0], 
									     amino_acids[index+2][1], 
									     amino_acids[index+2][2]+1) < 0) || 
	    (find_amino_acid_in_this_coordinate(amino_acids[index-1][0], 
						amino_acids[index-1][1], 
						amino_acids[index-1][2]-1) < 0 &&
	     find_amino_acid_in_this_coordinate(amino_acids[index+2][0], 
						amino_acids[index+2][1], 
						amino_acids[index+2][2]-1) < 0);
      }

      // if the crank position is not opposite, then turn is of course possible
      else
	okay_to_copy = true;

      // now we finally copy coords
      if (okay_to_copy) {
	// copy coords (index)
	available_next_coords[2*num_valid_crankshafts_possible][0] = crankshaft_possible_pairs[2*i][0];
	available_next_coords[2*num_valid_crankshafts_possible][1] = crankshaft_possible_pairs[2*i][1];
	available_next_coords[2*num_valid_crankshafts_possible][2] = crankshaft_possible_pairs[2*i][2];
      
	// copy coords (index+1)
	available_next_coords[(2*num_valid_crankshafts_possible)+1][0] = crankshaft_possible_pairs[(2*i)+1][0];
	available_next_coords[(2*num_valid_crankshafts_possible)+1][1] = crankshaft_possible_pairs[(2*i)+1][1];
	available_next_coords[(2*num_valid_crankshafts_possible)+1][2] = crankshaft_possible_pairs[(2*i)+1][2];   
      
	// increment valid pairs
	num_valid_crankshafts_possible++;
      }
    }
  }
  return num_valid_crankshafts_possible;
}

/* 
   This subproc finds out the common axis coords that (index-2), (index-1), (index+2), (index+3) have
   this subproc is to make crankshafts_possible look simpler
   THE ALGORITHM: if (index-1) and (index+2) are on the same Z coordinate (BUT (index) and (index+1)
   are not on that same Z coordinate), we try out combinations of 
   their neighbor coords for (index) and (index+1) as possible pairs for new crankshaft positions
   Repeat for X and Y axes
   ALTHOUGH COMPLEX, THIS CODE HAS BEEN TESTED 
*/
int crankshafts_possible_subproc(int index) {
  int code = 0, num_trial_crankshafts_possible = 0;
  if (amino_acids[index-1][0]  == amino_acids[index+2][0])
    code++;
  if (amino_acids[index-1][1]  == amino_acids[index+2][1])
    code += 3;
  if (amino_acids[index-1][2]  == amino_acids[index+2][2])
    code += 5;

  // make sure the proposed "crank" is not linear
  if (code == 4 && amino_acids[index-1][0]  == amino_acids[index][0] &&
      amino_acids[index-1][1]  == amino_acids[index][1])
    return 0;
  if (code == 6 && amino_acids[index-1][0]  == amino_acids[index][0] &&
      amino_acids[index-1][2]  == amino_acids[index][2])
    return 0;
  if (code == 8 && amino_acids[index-1][1]  == amino_acids[index][1] &&
      amino_acids[index-1][2]  == amino_acids[index][2])
    return 0;

  int i = 0;
  if (code == 4 || code == 6) {
    num_trial_crankshafts_possible += 2;
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
    num_trial_crankshafts_possible += 2;
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
    num_trial_crankshafts_possible += 2;
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
  /* for debugging purposes
  cout<<"possibles - "<< endl;
  for (int i=0; i<4; i++)
    cout << "(" << crankshaft_possible_pairs[2*i][0] << ", "
	 << crankshaft_possible_pairs[2*i][1] << ", "
	 << crankshaft_possible_pairs[2*i][2] << "), " 
	 << "(" << crankshaft_possible_pairs[(2*i)+1][0] << ", "
	 << crankshaft_possible_pairs[(2*i)+1][1] << ", "
	 << crankshaft_possible_pairs[(2*i)+1][2] << ")" << endl;
  cout<<endl;
  */
  return num_trial_crankshafts_possible;
}

// returns the number of coords possible for end amino acid to move to and 
// places the new coords into available_next_coords
int end_moves_possible(int index) {
  int num_end_moves_possible = 0;

  // in case function is called on a non-end protein, returns 0
  if (index != 0 && index != NUM_AMINO_ACIDS-1)
    return 0;

  // get neighbors
  if (index == 0)
    store_neighbors(1);
  else
    store_neighbors(NUM_AMINO_ACIDS-2);

  for (int j=0; j<6; j++) {
    // if there is no amino acid present in this coordinate, then coord into available_next_coords
    if (find_amino_acid_in_this_coordinate(neighbors[j][0], neighbors[j][1], neighbors[j][2]) < 0) {
      available_next_coords[num_end_moves_possible][0] = neighbors[j][0];
      available_next_coords[num_end_moves_possible][1] = neighbors[j][1];
      available_next_coords[num_end_moves_possible][2] = neighbors[j][2];
      num_end_moves_possible++;
    }
  }

  return num_end_moves_possible;
}

// save old coords into current_coords_being_displaced, then flip corner
void flip_corner(int index) {
  current_coords_being_displaced[0][0] = amino_acids[index][0];
  current_coords_being_displaced[0][1] = amino_acids[index][1];
  current_coords_being_displaced[0][2] = amino_acids[index][2];

  set_amino_acid_coords(index, available_next_coords[0][0], 
			available_next_coords[0][1], available_next_coords[0][2]);
  return;
}

// choose a crankshaft, save old coords into current_coords_being_displaced, then perform crankshaft
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

// save old coords into current_coords_being_displaced, then make end move
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
