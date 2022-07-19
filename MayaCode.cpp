//made in collboration with Dr.Paul Higgs team at McMaster Univeristy (c)2018-2019; in collboration with Vismay Shah and Andrew Tupper
// Globals
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>               
#include <vector>
#include <random>                
#include <stdlib.h>
#include <algorithm>
#include <cassert>

class Strand {
public:
	// We’ll store a few of the parameters as static as they are the same for all strands in the simulation
	static const std::string master_sequence;
	static const std::string master_structure;
	static double dk;
	static double dm;
	static double Mmin;
	static double Kmax;
	std::string sequence;
	double M;
	double K;
	double K_out;
	double M_out;
	double hs;
	double hd;
	std::string templates;
	Strand(FILE *ptr)
	{
	// constructor to initialize a strand using a FILE pointer -- for loading/resuming a simulation
	}

	Strand(std::string seq)
	{ // One of 2 regular constructors - only requires a sequence
	sequence = seq;
	M = calculate_M();
	K = calculate_K();
	}
	  
	double calculate_M()
	{
		//for ( int i=0; i< sequence.size(); i++){
			//for ( int j=0, j<mastersequence.size(); j++){
				//if( char sequence[i]== char master_sequence[j]){ //checks to see if the master template base is equal to the master sequence 
	/*int nucleo={10,20,30,40,50};	
	for (int i=0; i<nucleo.size(); i++){
		if (char sequence.nucleo[i]== char master_sequence.nucleo[i]){
			continue;
		}
		else{(hd=hd +1);		
			
		}  */
			//std::vector <Strand> nucleo = {10,20, 30,40,50};
			//if (char sequence.nucleo[i]== char master_sequence.nucleo[i])
		// Determine distance from this sequence to master sequence
	std::vector <int> nucleo(5);
	nucleo[0]=10;
	nucleo[1]=20;
	nucleo[2]=30;
	nucleo[3]=40;
	nucleo[4]=50;
		for (int i=0; i<nucleo.size(); i++){
			char baset = sequence[nucleo[i]];
			char basem = master_sequence[nucleo[i]];
			  if (baset==basem){
				    continue;
			  }
			  else{ (hd=hd +1);
			  }
		}
		return calculate_M(hd);
	}
	double calculate_M(double hd)
	{
	M = Mmin + dm * hd;
	return M;
	}
	
double calculate_K()
	{
	/*	// Determine distance from the sequence to compatible structure
			//for (int i=0; i<
		char bc='(';
		char be=')';
		char d= '.';
		std:: vector<char> bases;
		for( int i=0; i<sequence.size(); i++){
			char c1= sequence[i];
			char c2= sequence[-1*i];
				if (c1 == 'A') {
					if (c2 == 'U') {
						bases[i]=bc;
						bases[-1*i]=be;
				}
					}
				else if (c1 == 'U') {
					if (c2 == 'A') {
						bases[i]=bc;
						bases[-1*i]=be;
				}
					}
				else if (c1 == 'G') {
					if (c2 == 'C') {
						bases[i]=bc;
						bases[-1*i]=be;
				}
					}	 
				else if (c1 == 'C') {
					if (c2 == 'G') {
						bases[i]=bc;
						bases[-1*i]=be;
				}
					} 
				else { (bases[i]=d, bases[-1*i]=d);
		
				}
		}
	for (int j=0; j<bases.size(); j++){  //iterate through the bases vector 
			char seq_base = bases[j];
			
			char mas_base = master_structure[j];
			  if (seq_base==mas_base){
				    continue;
			   }
			   else{ ( hs=hs+1);
			   }
		  }	    */
	
	std::string baso= "............"
		for (int i=0; i<sequence.size(); i++){
				if 
				
				
	std::vector <int> hamm_struct=(6);
	hamm_struct[0]=3; //A
	hamm_struct[1]= 5; //A
	hamm_struct[2]= 4; //A
	hamm_struct[3]= 9; //U
	hamm_struct[4]= 11; //U
	hamm_struct[5]= 10; //U
	char bc='(';
		char be=')';
		char d= '.';
		for (int i=0; i<hamm_struct.size(); i++){
			char base_seq = sequence[hamm_struct[i]];
			char base_mast = master_sequence[hamm_struct[i]];
			  if (base_seq==base_mast){
				    if (i<2){
						baso[i]=bc
					};
					else {(baso[i]=be);
					}
			  }
			  else{ (continue);
			  }
		}
			
		return calculate_K(hs);
	}
	double calculate_K(double hs)
	{
		K= Kmax - dk * hs;
		return K;
	}
	void printStrand(FILE *ptr)
	{
			// Write out one strand to the provided output file
			// this should look something like:
		fwrite(sequence.c_str(), sizeof(char), sequence.size(), ptr);
		// ptr=fopen ("myfile.bin" , "wb");
		// fwrite(sequence.c_str(), sizeof(char), sizeof(sequence), ptr);
		// fclose(ptr);
		// but I haven’t tested it so try it out and try to fix it!
	}
};
const std::string Strand::master_sequence= "AGCAAAAUUUG";
const std::string Strand::master_structure= "..(((...)))";
double Strand::dk=4;
double Strand::dm=2;
double Strand::Mmin=4;
double Strand::Kmax=1;


class Site{
// Sites can hold multiple strands
// Replication stops in a site if S > Smax (a pre-defined constant)
public:
	static int Smax;	
	std::vector< Strand >  strands_in_site;  
	// Now we need constructors
	// First an empty one to just set up a lattice easily
	Site() {}
	// And one which can fully initialize the site
	Site( std::vector< Strand >  s);
	
	void addStrand(Strand t)
	{
		strands_in_site.push_back(t);
	}

	void printSite(FILE *ptr)
	{	
		int q= strands_in_site.size();
		// Write the vector size to the file
		 // ptr=fopen ("myfile.bin" , "wb");
		fwrite(&q, sizeof(q),1, ptr);
		// Iterate through the vector and call print on each strand
		// for s in strand (
		//	s.printStrand(ptr)
		for (int i=0; i<strands_in_site.size(); i++){
			 strands_in_site[i].printStrand(ptr);
		}
		// fclose(ptr);
	}

	Site(FILE* ptr)
	{
		// read in vector size
		unsigned int vect_size;
		fread(&vect_size, sizeof(unsigned int), 1, ptr);

		//ifstream seq_file(ptr)----I don't know if this is correct...
		// ptr=fopen("myfile.bin", "rb");
		// fread(s, sizeof(char), sizeof(s), ptr);
		// Iterate through vector size
		for (int i=0; i<vect_size; ++i) {
			Strand temp(ptr);
			strands_in_site.push_back(temp);
		}
		// for i in range vector_size
		//	Strand temp(ptr)
		//	strands_in_site.push_back(temp)
		// for (int i=0; i<len(s); i++){
		// 	Strand temp(ptr);
		// 	strands_in_site.pushback(temp);
			
		// }
		// fclose(ptr);
	}
};

int Site::Smax = 10;


// ============= ANDREWS CODE =====================

void deathStepSite(Site& site, float death_rate, float dt);
void replicateStepSite(Site& site, float dt);
Strand replicateStepStrand(const Strand& polymerase, const Strand& template0);
// Death functions
void printStats(const std::vector<std::vector<Site>>& lattice);

void deathStep(std::vector<std::vector<Site>>& lattice, float death_rate, float dt)
{ 
death_rate=1;
	// Iterate through every site and call deathStepSite
	// for i in len(lattice)
	for ( int i=0; i<lattice.size(); i++){
	// for ( int i=0; i<len(lattice[0]); i++){
	// 	for j in len(lattice[0])
		for ( int j=0; j<lattice[i].size(); j++){
	//		deathStepSite( lattice[i][j], ... )
			
		  deathStepSite(lattice[i][j], death_rate, dt);

		}
	}

}

void deathStepSite(Site& site, float death_rate, float dt)
{
  death_rate=1;
	// Create a new vector of type Strand called new_vect
	std::vector<Strand> new_vect;

	// Iterate through each strand in this site and call deathStepStrand on the strand
	for (int i=0; i<site.strands_in_site.size(); i++) {		
		// 	Death occurs with probability death_rate * dt which must be < 1
		// 	Determine random floating point number 'rf' in range [0,1)
		float rf = float(rand()) / RAND_MAX;

		// 	If rf < death_rate * dt
		if (rf < death_rate * dt) {
			// Strand is dead so don't add to new vector
			continue;
		}
		else {
			// Strand is still alive so add to new vector
			new_vect.push_back(site.strands_in_site[i]);
		}
	}

	// Set site.strands in site = new_vect
	site.strands_in_site = new_vect;
}

// Replicate functions
void replicateStep(std::vector<std::vector<Site>>& lattice, float dt)
{
	// Iterate through every site and call replicateStepSite
	// for i in len(lattice)
	for (int i=0; i<lattice.size(); i++) {
	// 	for j in len(lattice[i])
		for (int j=0; j<lattice[i].size(); j++) {	
		//		replicateStepSite( lattice[i][j], ... )
				replicateStepSite(lattice[i][j], dt);
		}
	}
}		

void replicateStepSite(Site& site, float dt)
{
	// Check if site has less strands then Smax
	if (site.strands_in_site.size() < site.Smax) {
		//	Get number of strands on this site num_strands
		int num_strands = site.strands_in_site.size();
		// 	Iterate through each strand in this site (for i in range(num_strands))
		for (int i=0; i<num_strands; i++) {
			//	Check if the strand has a non-zero rate of replication
			Strand polymerase = site.strands_in_site[i];
			if (polymerase.K != 0) {
				// Iterate through every other strand (excluding this one)
				for (int j=0; j<num_strands; j++) {
					Strand template0 = site.strands_in_site[j];
		
					if (i != j){ 
						  // Determine prob of replication
						  float prob_replication = polymerase.K * dt;
						  // Determine random floating point number 'rf' in range [0,1)
						  float rf = float(rand()) / RAND_MAX;
					
						  // If rf < prob_replication
						  if (rf < prob_replication){
						      site.strands_in_site.push_back(replicateStepStrand(polymerase,template0));
						  } //{-- for loading/resuming a simulation
					}
				
				}
			}
		}
	}
}

Strand replicateStepStrand(const Strand& polymerase, const Strand& template0)
{
	// Something is wrong here in terms of biology, what is it?
	// Determine m from polymerase strand
	float m = polymerase.M;

	// Create a blank string object
	// std::string new_seq
	std::string new_seq;

	// Iterate through each base in temp
	// for base in template
	// 	Check each case of base = A, U, G, C
	//	if base == A:
	//		generate random float 'rf' in range [0,1)d
	//		if rf < m: // mapped to A
	//			new_seq.push_back('A')
	//		else if rf < m + (1 - 3m):	// mapped to U
	//			new_seq.push_back('U')
	// 		else if rf < m + (1 - 3m) + m: // mapped to G
	//			new_seq.push_back('G')
	//		else:	// mapped to C
	//			new_seq.push_back('C')
	//	else if base == 'U':
	//		.
	//		.
	//		.
	// Done for loop
	// Create new Strand object using sequence and return
	// return Strand(new_seq)

	for (int i=0; i<template0.sequence.size(); i++){
		char base =template0.sequence[i];
		int m;
		if (base == 'A'){
			// generate random float 'rf' in range [0,1)d
			double rf = double(rand()) / RAND_MAX;
			if ( (0 <= rf) && (rf < m) ) { // mapped to A
				new_seq.push_back('A');
			}
			else if ( (m <= rf) && (rf < m + (1 - 3*m)) ){	// mapped to U
				new_seq.push_back('U');
			}
			else if ( (m + (1 - 3*m) <= rf) && (rf < m + (1 - 3*m) + m) ) { // mapped to G
				new_seq.push_back('G');
			}
			else if ( (m + (1 - 3*m) + m <= rf) && (rf < 1) ) { // mapped to C
				new_seq.push_back('C');
			}
			else { 
				assert(false);
			}
		}	
		
		else if (base == 'U'){
			char base = template0.sequence[i];
			for (int i=0; i<template0.sequence.size(); i++){
			// generate random float 'rf' in range [0,1)d
			double rf = double(rand()) / RAND_MAX;
			    if ( (0 <= rf) && (rf < m) ) { // mapped to A
				new_seq.push_back('U');
			    }
			    else if ( (m <= rf) && (rf < m + (1 - 3*m)) ){	// mapped to U
				new_seq.push_back('A');
			    }
			    else if ( (m + (1 - 3*m) <= rf) && (rf < m + (1 - 3*m) + m) ) { // mapped to G
				new_seq.push_back('G');
			    }
			    else if ( (m + (1 - 3*m) + m <= rf) && (rf < 1) ) { // mapped to C
				new_seq.push_back('C');
			    }
			    else { 
				assert(false);
			    }
			}
		}
		
		else if (base == 'C'){
			char base = template0.sequence[i];
			for (int i=0; i<template0.sequence.size(); i++){
			// generate random float 'rf' in range [0,1)d
			double rf = double(rand()) / RAND_MAX;
			    if ( (0 <= rf) && (rf < m) ) { // mapped to A
				new_seq.push_back('U');
			    }
			    else if ( (m <= rf) && (rf < m + (1 - 3*m)) ){	// mapped to U
				new_seq.push_back('G');
			    }
			    else if ( (m + (1 - 3*m) <= rf) && (rf < m + (1 - 3*m) + m) ) { // mapped to G
				new_seq.push_back('A');
			    }
			    else if ( (m + (1 - 3*m) + m <= rf) && (rf < 1) ) { // mapped to C
				new_seq.push_back('C');
			    }
			    else { 
				assert(false);
			    }
			}
		}
		
		else if (base == 'G'){
			char base = template0.sequence[i];
			for (int i=0; i<template0.sequence.size(); i++){
			// generate random float 'rf' in range [0,1)d
			double rf = double(rand()) / RAND_MAX;
			    if ( (0 <= rf) && (rf < m) ) { // mapped to A
				new_seq.push_back('U');
			    }
			    else if ( (m <= rf) && (rf < m + (1 - 3*m)) ){	// mapped to U
				new_seq.push_back('C');
			    }
			    else if ( (m + (1 - 3*m) <= rf) && (rf < m + (1 - 3*m) + m) ) { // mapped to G
				new_seq.push_back('A');
			    }
			    else if ( (m + (1 - 3*m) + m <= rf) && (rf < 1) ) { // mapped to C
				new_seq.push_back('G');
			    }
			    else { 
				assert(false);
			    }
			}
		}
	}
	// Reverse the string
	//std::string new_new_seq = new_seq.reverse();
	std::reverse(new_seq.begin(), new_seq.end());

	// Create new Strand instance and return 
	return Strand(new_seq);
	

}

// Getting stats for printing
void printStats(const std::vector<std::vector<Site>>& lattice)
{
	// Variables to update:
	float m_sum = 0;
	float m_sqr_sum = 0;
	float k_sum = 0;
	float k_sqr_sum = 0;
	float num_parasites=0 ;
	float num_strands= 0;
	// Iterate through every site in lattice
	//for i in len(lattice)
	//for j in len(lattice[i]))
	for (int i=0; i<lattice.size(); i++){
		for(int j=0; j<lattice[i].size(); j++){
		    for( int k=0; k<lattice[i][j].strands_in_site.size(); k++){
		      Strand sequence_strandk = lattice[i][j].strands_in_site[k];
		      float k2 = sequence_strandk.K;
		      float m2 = sequence_strandk.M; 
		      m_sum= m_sum + m2;
		      k_sum= k_sum + k2;
		
		      
		    }
		}
	}	
	for (int i=0; i<lattice.size(); i++){
	    for(int j=0; j<lattice.size(); j++){
		    int num_strands0 = lattice[i][j].strands_in_site.size();
		    float num_strands= num_strands + num_strands0;
	    }
	}
	
	for (int i=0; i<lattice.size(); i++){
	    for( int j=0; j<lattice[i].size(); j++){
			for (int k=0; k< lattice[i][j].strands_in_site.size(); k++){
				Strand num_parasites0= lattice[i][j].strands_in_site[k];
					if ( num_parasites0.K = 0 ){
						num_parasites = num_parasites + 1;
					}
			}			
	    }
	}
	 
	//Get k and m of strand
	//		add to sums and sqr sums
	 float m_avg = m_sum / num_strands;
	float m_var = m_sqr_sum - (m_avg * m_avg);
	// Do same for k
	float k_avg= k_sum / num_strands;
	float k_var=k_sqr_sum - (k_avg*k_avg);
	std::cout<< "Mutation average:"<< m_avg<< ','<< "Mutation variance:" <<m_var<< ','<< "Replication average:"<< k_avg<< ','<<"Replication variance"<< k_var << ','<<"NUmber of Strands:"<< num_strands<<std::endl;
	// Then print to these stats to the screen
};

// ==============================================



int main(int argc, char* argv[]) { // the additional params are for command line inputs
 
	
	// Start by defining some key parameters
	int length = 50;// some number
	int num_strands_init = 3;// some number - this is the number of polymerase sequences to be placed \
	in each Site upon initialization
	float time_max = 1;// some value
	float dt = 0.001;// some reasonably small value
	FILE *out; // pointer to the output file
	std::string outname = "default_output_file.txt";
	float death_rate=1;
	int l_size=50;
	// If we’re using command line args, replace the default values with the ones provided
	if (argc > 1){
		length = std::atoi(argv[1]);
		num_strands_init = std::atoi(argv[2]);
		time_max = std::atof(argv[3]);
		dt = std::atof(argv[4]);
		outname = argv[5];
	}
	out = fopen(outname.c_str(), "w"); // open the output file in write mode

	int sim_t_max = int(time_max / dt);

	// Now we want to set up a blank lattice, and generate a Strand object from a master sequence
	std::vector< std::vector <Site> > lattice; // initialize blank lattice	
	std::string poly_seq = "AGGGGUCCCUA"; // fill this with a string from Andrew’s file
	//std::vector< Strand > fill_site;

	// iterate from 0 to num_strands_init and allocate that many polymerases to each site
	for( int i=0; i<l_size; i++){
		lattice.push_back(std::vector<Site>());
		for(int j=0; j<l_size; j++){
			lattice[i].push_back(Site());
			for (int p = 0; p < num_strands_init; p++){
				//for( int g=0; g<poly_seq.size(); g++){
				// We start out by creating ‘perfect’ polymerases, with hd=0 and hs=0
				Strand polymerase(poly_seq);
				//if (polymerase.hd !=0 && polymerase.hs !=0){
				// fill_site.push_back(     //create the strands here, make sure to use the right constructor!    );
				//lattice[i][j].strands_in_site.push_back(polymerase);		  
				lattice[i][j].addStrand(polymerase);
				std::cout<< latticeStrands<<lattice<<std::endl; 
			}
		}	
	}

	// Set up for a way to easily visit Sites in a random order
	std::vector< std::pair<int, int> > coords;
	// And now we properly set up the Sites by filling them with the Strands
	// iterate over each site
	// for (int i …{
	for (int i=0; i< lattice.size();i++){
	// for (int j…{
		for (int j=0; j< lattice[i].size();j++){

			// and fill the site aka call the appropriate site constructor	
			// 	lattice[i][j] = …
			//lattice[i][j]= .strands_in_site; 
			// add the current lattice coordinate to the list of coordinates
			coords.push_back(std::make_pair(i,j)); 

		}
	}
	
	
				
		
	
		// Main program loop
	for (int sim_t = 0; sim_t <sim_t_max; sim_t++){
		// Shuffle the coordinates vector to get a random visitation order for the Sites
		replicateStep(lattice, dt);
		deathStep(lattice,death_rate,dt);
		std::random_shuffle(coords.begin(), coords.end());
		for (int i = 0; i < coords.size(); i++){
			// Get the x,y of the lattice site
			int x=coords[i].first; 
			int y=coords[i].second; 
			// Call replication, death and diffusion here for site x,y
			//lattice[x][y].No diffusion as of yet"
			
		}
		printStats(lattice);
	};
// Deal with the statistical calculations here

	fclose(out); // close the file
	return 0;
}
