#include <vector>
#include <random>
#include "time.h"
#include <iostream>
#include <stdio.h>
#include <fstream>


//create a sequence of 33 basepairs such that all 16 bp steps included twice. The order of the steps is random. The seed of the random generator is the clock time.


using namespace std;

int A = 0;
int T = 1;
int C = 2;
int G = 3;

int N;
int seed;

vector<int> bases;
ofstream ofile;

int main(int argc, char **argv) {

cout << "Generating random (random bases) double helix" << endl;

if(argc!=3) {
	cout << "Something is not ok" << endl;
	cout << "Usage: " << endl;
	cout << argv[0] << " Number_of_bases seed" << endl;
	return 1;
}

N = atoi(argv[1]);
seed = atoi(argv[2]);


/*setup rng mt19937*/
mt19937 rng(seed);	//seed = time
uniform_real_distribution<double> udist01(0,1);

//create a sequence

for(int i = 0; i < N; i++) {
	bases.push_back((int) 4*udist01(rng));
}

//cout << "Generated!" << endl;

char type;
ofile.open("seq.txt");
ofile << "DOUBLE ";
for(int i = 0; i < N; i++) {
	if(bases[i] == A) type = 'A';
	else if(bases[i] == T) type = 'T'; 
	else if(bases[i] == C) type = 'C'; 
	else if(bases[i] == G) type = 'G'; 
	ofile << type;
}

ofile.close();

return 0;

}
