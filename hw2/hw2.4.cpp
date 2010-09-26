#include "common.h"
#include "common.cpp"

#define MY_SEED 279683937

// declare the Point datatype - holds the coordinate and the array index of the bead that lives in it
typedef struct {
  int x, y, z, index;
} Point;

// declare functions for use
int find_bead_in_this_coordinate(int x, int y, int z);
double energy();

// global variables - very bad programming style, but code more optimized
const int NUM_AMINO_ACIDS = 36;
const char AMINO_ACID_TYPE[NUM_AMINO_ACIDS+1] = "HHPHHPPHPPPHPPHPHHHHPPPPHHPPHPHHHHPP";
int bead_positions[NUM_AMINO_ACIDS][3];
vector<Point> coordinates_containing_beads;


int main( int argc, char **argv ) {
  vector<int> a;
  a.push_back(2);
  a.push_back(5);
  vector<int>::iterator it = find(a.begin(),a.end(),7);
  //if (  it==0)
  //cout<<"true"<<endl;
  //  else
    cout<< "false" <<endl;
}


double energy() {
  double energy = 0.0;
  int index, neighbors[6][3];
  for (int i=0; i<NUM_AMINO_ACIDS; i++) {
    // stores coordinates of immediate neighbors in array
    neighbors[0][0] = bead_positions[i][0]+1;
    neighbors[0][1] = bead_positions[i][1];
    neighbors[0][2] = bead_positions[i][2];

    neighbors[1][0] = bead_positions[i][0]-1;
    neighbors[1][1] = bead_positions[i][1];
    neighbors[1][2] = bead_positions[i][2];

    neighbors[2][0] = bead_positions[i][0];
    neighbors[2][1] = bead_positions[i][1]+1;
    neighbors[2][2] = bead_positions[i][2];

    neighbors[3][0] = bead_positions[i][0];
    neighbors[3][1] = bead_positions[i][1]-1;
    neighbors[3][2] = bead_positions[i][2];

    neighbors[4][0] = bead_positions[i][0];
    neighbors[4][1] = bead_positions[i][1];
    neighbors[4][2] = bead_positions[i][2]+1;

    neighbors[5][0] = bead_positions[i][0];
    neighbors[5][1] = bead_positions[i][1];
    neighbors[5][2] = bead_positions[i][2]-1;

    // for each of the neigbor positions, ind the index of bead that lives in it, if any, and add inter-bead energy to total energy
    for (int j=0; j<6; j++) {
      index = find_bead_in_this_coordinate(neighbors[i][0], neighbors[i][1], neighbors[i][2]);
      
      if (index >= 0) {
	if (AMINO_ACID_TYPE[i] == 'H')
	  energy += (AMINO_ACID_TYPE[index] == 'H') ? -1.0 : 0.0;
	else
	  energy += (AMINO_ACID_TYPE[index] == 'H') ? 0.0 : -0.5;
      }

      else
	energy += (AMINO_ACID_TYPE[i] == 'H') ? 0.5 : -0.5;
    }
  }
  
  // accounts for double-counting
  return energy / 2.0;
}

// finds the array index of the bead that lives in coordinates (x,y,z) or returns -1 otherwise
int find_bead_in_this_coordinate(int x, int y, int z) {
  for (vector<Point>::iterator i=coordinates_containing_beads.begin(); 
       i!=coordinates_containing_beads.end();
       i++) {
    if ((*i).x == x && (*i).y == y && (*i).z == z)
      return (*i).index;
  }
  return -1;
}
