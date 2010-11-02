/*
	Monte Carlo simulation of lattice protein
	Original by: Enghui Yap, Feb 08
	Modified by: Aaron Kaluszka, Feb 09, Shachi Katira, Feb '10

	To compile, use the following command:
		g++ -lrt -O3 -o mclattice mclattice_student.cpp
	To run the program, use the following command:
		./mclattice <parameter filename> <config filename>
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <cmath>
#include <string.h>
#include "point.h"
#include "urng.h"

using namespace std;

typedef vector<Point> PointV;
#define EMPTY 255

#define SOLVENT 1 //SET SOLVENT 1 FOR INTERACTION MATRICES WITH HPS

void readParams(char* pfname, double& temp, int& ncycles, int& ntraj, int& nSampleFreq, int& neqstep);
void readConfig(char* cfname, string& beadType, PointV& c, double V[]);
double energy(const string& beadType, const PointV& c, double V[]);
double updateEnergy(const string& beadType, const PointV& c, double V[], double potential);
bool mcmoves(const PointV& cn, int jm);
bool mcaccept(double vcurr, double vtrial, double beta);
void print_pos(const PointV& c, string& beadType);
bool move_endbeads(int j, const PointV& c);
bool move_corner(int j, const PointV& c);
bool move_crank(int j, const PointV& c);
int idum;
unsigned char space[16777216]; // lookup table
Point changes[2];
int changesindex[2];

// ==================================================
// main program 
// ==================================================
int main(int argc, char** argv )
{
	// monte carlo variables
	double vcurr, vtrial, vmin;
	double beta;

	// simulation parameters
	int ncycles; // no. of monte carlo steps to run for each trajectory
	int ntraj; // no. of independent trajectories to launch for statistics
	int nSampleFreq; // no. of steps between sampling points
	double kT; // temperature in reduced units

	// beadtype and config variables
	string beadType; 
	PointV ci, cr;
	int nbeads;
	double V[65536];

	// variables for averages
	int nsample;
	int neqstep; // equilibration steps before evaluating averages
	int bead;
	double vsum, vsqsum, v_avg, vsq_avg, cv;


	//timespec start, end;
	double diff;

	char* pfname = argv[1]; // input argument 1: parameter filename
	char* cfname = argv[2]; // input argument 2: config filename

	idum = time(NULL); // initialize random seed

	// read in simulation parameters
	readParams(pfname, kT, ncycles, ntraj, nSampleFreq, neqstep);

	// read in beadtype and initial configuration
	readConfig(cfname, beadType, ci, V);
	nbeads = ci.size();

	//===================================================
	// ACTUAL MONTE CARLO RUN STARTS HERE

	beta = 1.0 / kT;
	ofstream avgfile("averagesolv.out"); // for thermodynamic averages

        clock_t start = clock();

	//clock_gettime(CLOCK_MONOTONIC, &start);

	for(int r = 0; r < ntraj; r++)
	{
		vsum = vsqsum = v_avg = vsq_avg = cv = nsample = 0; // initialize averages to zero
		cr = ci;// initialize config for this trajectory to ci
		memset(space, EMPTY, 16777216);
		for(int i = 0; i < nbeads; i++) space[ci[i].get()] = i;

		//evaluate potential of starting configuration
		vmin = vcurr = energy(beadType, cr, V);
		cout << "Starting Potential for Monte Carlo: " << vcurr << endl;

		for(int n = 0; n < ncycles; n++)
		{
			for(int m = 0; m < nbeads; m++)
			{
				changesindex[0] = changesindex[1] = -1;
				if(mcmoves(cr, m))
				{
					vtrial = updateEnergy(beadType, cr, V, vcurr);
					if(mcaccept(vcurr, vtrial, beta))
					{
							
						vcurr = vtrial;
						for(int i = 0; i < 2 && changesindex[i] > -1; i++)
						{
							
							bead = changesindex[i];
							space[cr[bead].get()] = EMPTY;
							space[changes[i].get()] = bead;
							cr[bead] = changes[i];
						}
						
						
					}
				}
			}
			if(vtrial < vmin)  vmin = vcurr;

			if(n % nSampleFreq == 0 && n >= neqstep) 
			{
				nsample++; // no.of sample points

				//========== START STUDENT INPUT ==========
				// Calculate the running statistical averages 
				// You can use any of these variables to save values: 
				//	vsum, vsqsum, v_avg, vsq_avg, cv
				// You can use these variables to do your calculations:
				//	vcurr, beta (which is 1/kT)
				
				//=========== END STUDENT INPUT ===========

				// writes out: no. of samples, <V>, <V^2>, Cv to avgfile
				//avgfile << nsample << "\t" << v_avg << "\t" << vsq_avg << "\t" << cv << endl;
			}
		}
		cout << "Minimum Potential from Monte Carlo: " << vmin << endl;
		
		// writes out: no. of samples, <V>, <V^2>, Cv to avgfile
		avgfile << r << "\t" << v_avg << "\t" << vsq_avg << "\t" << cv << endl;
	}

        clock_t end = clock();

	//clock_gettime(CLOCK_MONOTONIC, &end);

	// Compute how long it takes to run the code

        diff = (double(end) - double(start))/CLOCKS_PER_SEC;
	//diff = (double(end.tv_sec - start.tv_sec) * 1000000000 + end.tv_nsec - start.tv_nsec) / 1000000000;
	cout << "Program run time: " << diff << " seconds" << endl;

  	avgfile.close();
	return 0;
}

// define functions here
void readParams(char* pfname, double& temp, int& ncycles, int& ntraj, int& nSampleFreq, int& neqstep)
{
	ifstream pfile(pfname);
	if(!pfile) 
	{
		cerr << "Error: cannot open " << pfname << endl;
		exit(1);
	}

	pfile >> temp // temp in reduced units
		>> ncycles // no. of steps to run at each temp
		>> ntraj // no. of trajectories to collect
		>> nSampleFreq // no. of steps between samples
		>> neqstep; // Equilibration steps before collecting statistics

	cout << "Simulation Parameters: " << endl
		<< " Temperature: " << temp << endl
		<< " No. of steps: " << ncycles << endl
		<< " No. of trajectories: " << ntraj << endl
		<< " Steps between samples: " << nSampleFreq << endl
		<< " Equilibration steps before collecting statistics: " << neqstep << endl;
}

// read in beadtypes and initial configuration
void readConfig(char* cfname, string& beadType, PointV& c, double V[])
{
	int x, y, z;
	char btype;
	string beadtypes;
	int ntypes;
	double energy;
	int nbeads;

	ifstream configfile(cfname);
	if(!configfile) 
	{
		cerr << "Error: cannot open " << cfname << endl;
		exit(1);
	}
	configfile >> beadtypes;
	ntypes = beadtypes.size();
	for(int i = 0; i < ntypes; i++)
	{
		for(int j = 0; j < ntypes; j++)
		{
			configfile >> energy;
			V[((unsigned int)beadtypes[i] << 8) | beadtypes[j]] = energy;
		
		
		}

	}
	

	configfile >> nbeads;

	for(int i = 0; i < nbeads; i++)
	{
		configfile >> btype >> x >> y >> z;
		c.push_back(Point(x, y, z));
		beadType.push_back(btype);
	}

	// Prints out for debugging
	cout << "Read following chain configuration:" << endl;
	for(int i = 0; i < nbeads; i++)
		cout << "Bead #" << i << "\t" << beadType[i] << " " << c[i] << endl;
	configfile.close();
}

// returns the energy of the protein based on its config and bead sequence
double energy(const string& beadType, const PointV& c, double V[])
{
	// Loop through all unique, non-bonded i-j pairs: if they 
	// are neighbours then add up their energies
	unsigned char bead2;
	double potential = 0;
	int nbeads = c.size();
	static Point neighbors[6] = {Point(1, 0, 0), Point(-1, 0, 0), Point(0, 1, 0), 
		Point(0, -1, 0), Point(0, 0, 1), Point(0, 0, -1)};

	for(int i = 0; i < nbeads; i++)
	{
		for(int j = 0; j < 6; j++)
		{
			bead2 = space[(c[i] + neighbors[j]).get()];
			if(bead2 != EMPTY && abs(i - bead2) > 1)
				potential += V[((unsigned int)beadType[i] << 8) | beadType[bead2]];
#ifdef SOLVENT==1
			if(bead2 == EMPTY)	
				potential += 2*V[((unsigned int)beadType[i] << 8) | 'S']; 
#endif				
		}
	}
	return potential / 2;
}

double updateEnergy(const string& beadType, const PointV& c, double V[], double potential)
{
	int bead1, bead2, bead3;
	static Point neighbors[6] = {Point(1, 0, 0), Point(-1, 0, 0), Point(0, 1, 0), 
		Point(0, -1, 0), Point(0, 0, 1), Point(0, 0, -1)};

	for(int i = 0; i < 2 && changesindex[i] > -1; i++)
	
	{
		bead1 = changesindex[i];

		for(int j = 0; j < 6; j++)
		{
			bead2 = space[(c[bead1] + neighbors[j]).get()];
			bead3 = space[(changes[i] + neighbors[j]).get()];
	
			if(bead2 != EMPTY && abs(bead1 - bead2) > 1)
				potential -= V[((unsigned int)beadType[bead1] << 8) | beadType[bead2]]; //discount current neighbour bead interaction
			
			if(bead3 != EMPTY && abs(bead1 - bead3) > 1) 
				potential += V[((unsigned int)beadType[bead1] << 8) | beadType[bead3]]; //count new neighbour  bead interaction
	
#ifdef SOLVENT==1
			if(bead3 == EMPTY)
				potential += V[((unsigned int)beadType[bead1] << 8) | 'S']; 
			
			if(bead2 == EMPTY)
				potential -= V[((unsigned int)beadType[bead1] << 8) | 'S']; 
				
			if(bead2 != EMPTY && abs(bead1 - bead2) == 1)
				potential -= V[((unsigned int)beadType[bead2] << 8) | 'S'];
			
			if(bead3 != EMPTY && abs(bead1 - bead3) == 1)
				potential += V[((unsigned int)beadType[bead3] << 8) | 'S'];
				
			if(bead3 != EMPTY && abs(bead1 - bead3) > 1) 
				potential -= V[((unsigned int)beadType[bead3] << 8) | 'S'];  //discount new neighbour bead-solvent interaction
			
			if(bead2 != EMPTY && abs(bead1 - bead2) > 1)
				potential += V[((unsigned int)beadType[bead2] << 8) | 'S'];  //count incumbent neighbour bead - solvent interaction
#endif			
			
			
				
		 	
	    	}
	}
	return potential;
}

//MCMOVES : executes either Endbeads, Corner or Crankshaft move on a random bead
bool mcmoves(const PointV& cn, int jm)
{
	// Find the possible moveset
	// Faster than randomly selected moveset
	// This satisfies detailed balance because each bead
	// has only one possible movetype. So the new config
	// has the same available movetype as the old config
	return move_endbeads(jm, cn) || move_corner(jm, cn) || move_crank(jm, cn);
}

// MCACCEPT : 
// - computes the boltzman weight of the trial potential
// - compares to boltzman weight of current potential 
// - decides whether to accept
// c is the actual config and cn is the trial config
bool mcaccept(double vcurr, double vtrial, double beta)
{
	double vdiff, bzm, urn;
	bool accept = false;
	//========== START STUDENT INPUT ==========
	// Code your metropolis acceptance criteria here
	// You can use ran3(idum) for your urn
	// You can use vdiff to store delta-V and bzm to store the Boltzman factor
	// If move accepted, set accept to true
	// beta is equal to 1/kT


	//=========== END STUDENT INPUT ===========
	return accept;
}

// PRINT_POS : prints config to screen for checking
void print_pos(const PointV& c, string& beadType)
{ 
	int n = c.size();
	for(int i = 0; i < n; i++) cout << beadType[i] << " "  <<  c[i] << endl;
}

// MOVE_ENDBEADS: makes a move on the specified end bead
bool move_endbeads(int j, const PointV& c)
{
	int jend = 0, jend_nxt = 1;
	int endrand;
	int rEndmove[4];
	Point rEndmove_tmp[6];
	int nbeads = c.size();

  if(j == nbeads - 1)
	{
		jend = nbeads - 1;
		jend_nxt = nbeads - 2;
	}
	else if(j != 0) return false;

	// list all the 6 possible moves pivoting at jend_nxt
	rEndmove_tmp[0].set(c[jend_nxt].x - 1, c[jend_nxt].y,     c[jend_nxt].z    );
	rEndmove_tmp[1].set(c[jend_nxt].x + 1, c[jend_nxt].y,     c[jend_nxt].z    );
	rEndmove_tmp[2].set(c[jend_nxt].x,     c[jend_nxt].y - 1, c[jend_nxt].z    );
	rEndmove_tmp[3].set(c[jend_nxt].x,     c[jend_nxt].y + 1, c[jend_nxt].z    );
	rEndmove_tmp[4].set(c[jend_nxt].x,     c[jend_nxt].y  ,   c[jend_nxt].z - 1);
	rEndmove_tmp[5].set(c[jend_nxt].x,     c[jend_nxt].y  ,   c[jend_nxt].z + 1);

  // eliminate moves that clash with existing beads  
	int endmovetype = 0;
	for(int i = 0; i < 6; i++)
	{
		if(space[rEndmove_tmp[i].get()] == EMPTY)
		{
			rEndmove[endmovetype] = i;
			endmovetype++;
		}
	}
  
  // choose randomly from the possible end moves
  if(endmovetype > 0)
	{
		endrand = int(ran3(idum) * endmovetype);
		changesindex[0] = jend;
		changes[0] = rEndmove_tmp[rEndmove[endrand]];
		for(int exe=0; exe < endmovetype; exe++)
		return true;
	}
  return false;
}

//MOVE_CORNER : executes Corner move on a selected bead if possible
bool move_corner(int j, const PointV& c)
{
  int nbeads = c.size();
	if(j==0 || j == (nbeads - 1))
		return false;
// new position of corner bead
  int xc = c[j].x;
	int yc = c[j].y;
	int zc = c[j].z;
  
	// flips y and z 
	if(c[j + 1].x - c[j - 1].x == 0) 
	{
		if(c[j + 1].y - c[j].y == 0) 
		{
			yc = c[j - 1].y;
			zc = c[j + 1].z;
		}
		else 
		{
			yc = c[j + 1].y;
			zc = c[j - 1].z;
		}
	}
  
	// flips x and z
	else if(c[j + 1].y - c[j - 1].y == 0) 
	{
		if(c[j + 1].z - c[j].z == 0) 
		{
			zc = c[j - 1].z;
			xc = c[j + 1].x;
		}
		else 
		{
			zc = c[j + 1].z;
			xc = c[j - 1].x;
		}
	}
  
	// flips x and y
	else if(c[j + 1].z - c[j - 1].z == 0) 
	{
		if(c[j + 1].x - c[j].x == 0) 
		{
			xc = c[j - 1].x;
			yc = c[j + 1].y;
		}
		else 
		{
			xc = c[j + 1].x;
			yc = c[j - 1].y;
		}
	}
	else return false;
  
	Point pt(xc, yc, zc);
	if(space[pt.get()] == EMPTY)
	{
		changesindex[0] = j;
		changes[0] = pt;
		return true;
	}
	return false;
}

//MOVE_CRANK : executes Crankshaft move on specified bead
bool move_crank(int j, const PointV& c)
{
	char const_plane = 'N'; // initialize const plane to 'No'
	int janchor1, janchor2, jturn1, jturn2; // indices of anchor beads and turning beads
	Point rTurn1_tmp[4], rTurn2_tmp[4];
	int nbeads = c.size();
	if(j == 0 || j == (nbeads - 1))
		return false;
	

	Point dr = c[j + 2] - c[j - 1];
	if(dr.normsq() == 1) 
	{
		janchor1 = j - 1;
		janchor2 = j + 2;
		jturn1   = j;
		jturn2   = j + 1;
		if(dr.x != 0) const_plane = 'X';
		else if(dr.y != 0) const_plane = 'Y';
		else if(dr.z != 0) const_plane = 'Z';

	}
	else 
	{
		dr = c[j + 1] - c[j - 2];
		if(dr.normsq() == 1)
		{
			janchor1 = j - 2;
			janchor2 = j + 1;
			jturn1   = j - 1;
			jturn2   = j;	
			if(dr.x != 0) const_plane = 'X';
			else if(dr.y != 0) const_plane = 'Y';
			else if(dr.z != 0) const_plane = 'Z';
		}
	}
 	

	if(jturn1==0 || jturn2==0 || jturn1 == nbeads - 1 || jturn2 == nbeads -1)
		return false;

	switch(const_plane)
	{
	case 'X':
		//list all 4 rotations possible - rotate about anchors in Y-Z plane
		rTurn1_tmp[0].set(c[janchor1].x, c[janchor1].y - 1, c[janchor1].z    );
		rTurn2_tmp[0].set(c[janchor2].x, c[janchor2].y - 1, c[janchor2].z    );
		rTurn1_tmp[1].set(c[janchor1].x, c[janchor1].y + 1, c[janchor1].z    );
		rTurn2_tmp[1].set(c[janchor2].x, c[janchor2].y + 1, c[janchor2].z    );
		rTurn1_tmp[2].set(c[janchor1].x, c[janchor1].y,     c[janchor1].z - 1);
		rTurn2_tmp[2].set(c[janchor2].x, c[janchor2].y,     c[janchor2].z - 1);
		rTurn1_tmp[3].set(c[janchor1].x, c[janchor1].y,     c[janchor1].z + 1);
		rTurn2_tmp[3].set(c[janchor2].x, c[janchor2].y,     c[janchor2].z + 1);
		break;

	case 'Y':
		//list all 4 rotations possible - rotate about anchors in X,Z plane
		rTurn1_tmp[0].set(c[janchor1].x - 1, c[janchor1].y, c[janchor1].z    );
		rTurn2_tmp[0].set(c[janchor2].x - 1, c[janchor2].y, c[janchor2].z    );
		rTurn1_tmp[1].set(c[janchor1].x + 1, c[janchor1].y, c[janchor1].z    );
		rTurn2_tmp[1].set(c[janchor2].x + 1, c[janchor2].y, c[janchor2].z    );
		rTurn1_tmp[2].set(c[janchor1].x,     c[janchor1].y, c[janchor1].z - 1);
		rTurn2_tmp[2].set(c[janchor2].x,     c[janchor2].y, c[janchor2].z - 1);
		rTurn1_tmp[3].set(c[janchor1].x,     c[janchor1].y, c[janchor1].z + 1);
		rTurn2_tmp[3].set(c[janchor2].x,     c[janchor2].y, c[janchor2].z + 1);
		break;

	case 'Z':
		//list all 4 rotations possible - rotate about anchors in X,Y plane
		rTurn1_tmp[0].set(c[janchor1].x - 1, c[janchor1].y,     c[janchor1].z);
		rTurn2_tmp[0].set(c[janchor2].x - 1, c[janchor2].y,     c[janchor2].z);
		rTurn1_tmp[1].set(c[janchor1].x + 1, c[janchor1].y,     c[janchor1].z);
		rTurn2_tmp[1].set(c[janchor2].x + 1, c[janchor2].y,     c[janchor2].z);
		rTurn1_tmp[2].set(c[janchor1].x,     c[janchor1].y - 1, c[janchor1].z);
		rTurn2_tmp[2].set(c[janchor2].x,     c[janchor2].y - 1, c[janchor2].z);
		rTurn1_tmp[3].set(c[janchor1].x,     c[janchor1].y + 1, c[janchor1].z);
		rTurn2_tmp[3].set(c[janchor2].x,     c[janchor2].y + 1, c[janchor2].z);
		break;

	default: return false;
  }

	int crankmovetype, randcrank;
	int rTurn[3];
   
	//select only allowable moves that do no clash with existing beads
	crankmovetype = 0;
	for(int i = 0; i < 4; i++)
	{
		if(space[rTurn1_tmp[i].get()] == EMPTY && space[rTurn2_tmp[i].get()] == EMPTY)
		{
			rTurn[crankmovetype] = i;
			crankmovetype++;
		}
	}

	//choose a random allowable crankshaft move if there is any
	if(crankmovetype > 0)
	{
		randcrank = int(ran3(idum) * crankmovetype);
		changesindex[0] = jturn1;
		changes[0] = rTurn1_tmp[rTurn[randcrank]];
		changesindex[1] = jturn2;
		changes[1] = rTurn2_tmp[rTurn[randcrank]];
		return true;
	}
	return false;
}
