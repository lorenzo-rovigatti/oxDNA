#include <vector>
#include <random>
#include "time.h"
#include <iostream>
#include <fstream>
#include <stdio.h>


//create a sequence of 33 basepairs such that all 16 bp steps included twice. The order of the steps is random. The seed of the random generator is the clock time.


using namespace std;

int A = 0;
int T = 1;
int C = 2;
int G = 3;


struct BasePair {
	int ty1;
	int ty2;
};

vector<BasePair*> bps;
int* indices;
vector<BasePair*> ordered_bps;
vector<int> compatible_bps;
int id;
BasePair * bp;
ofstream ofile;

int seed;

int main(int argc, char**argv) {

seed = atoi(argv[1]);
cout << "Generating 32 bps (symmetry not accounted for) double helix with all steps twice." << endl;

//all possible junctions, twice.
for(int k = 0; k < 2; k++) {
	for(int i = 0; i < 4; i++) {
		for(int j = 0; j < 4; j++) {
			bp = new BasePair;
			bp->ty1 = i;
			bp->ty2 = j;
			bps.push_back(bp);
		}
	}
}

char type1,type2;
for(int i = 0; i < bps.size(); i++) {
	if(bps[i]->ty1 == A) type1 = 'A';
	else if(bps[i]->ty1 == T) type1 = 'T'; 
	else if(bps[i]->ty1 == C) type1 = 'C'; 
	else if(bps[i]->ty1 == G) type1 = 'G';
	if(bps[i]->ty2 == A) type2 = 'A';
	else if(bps[i]->ty2 == T) type2 = 'T'; 
	else if(bps[i]->ty2 == C) type2 = 'C'; 
	else if(bps[i]->ty2 == G) type2 = 'G'; 
	//cout << i << " " << type1 << " " << type2 << endl;
}

/*setup rng mt19937*/
mt19937 rng(seed);
uniform_real_distribution<double> udist01(0,1);

//create a sequence
int N = bps.size();
indices = new int[N];

int counts = 0;
while(ordered_bps.size()<bps.size()) {

	//if(counts == 1000) break;
	
	if(N == bps.size()) {
		for(int i = 0; i < bps.size(); i++) indices[i] = i;
		id = (int) N*udist01(rng);
		//cout << "id: " << id << endl;
		ordered_bps.push_back(bps[id]);
		for(int i = id; i < bps.size()-1; i++) {
			indices[i] = indices[i+1];
		}
		indices[bps.size()-1]=-1;
		for(int i = 0; i < bps.size(); i++) {
			//cout << " " << indices[i];
		}
		//cout << endl;
		N--;	
	}
	
	compatible_bps.clear();
	//cout << "obp ty2: " << ordered_bps[ordered_bps.size()-1]->ty2 << endl;
	for(int i = 0; i < N; i++) {
		if(bps[indices[i]]->ty1 == ordered_bps[ordered_bps.size()-1]->ty2) {
			//cout <<"Compatible: " << indices[i] << endl;
			compatible_bps.push_back(i);
		}
	}
	if(compatible_bps.size()==0) {
	
		//cout << "Attempt: " << counts << endl;
		//cout << "N: " << bps.size()-N << endl;
		for(int i = 0; i < ordered_bps.size(); i++) {
			if(ordered_bps[i]->ty1 == A) type1 = 'A';
			else if(ordered_bps[i]->ty1 == T) type1 = 'T'; 
			else if(ordered_bps[i]->ty1 == C) type1 = 'C'; 
			else if(ordered_bps[i]->ty1 == G) type1 = 'G';
			if(ordered_bps[i]->ty2 == A) type2 = 'A';
			else if(ordered_bps[i]->ty2 == T) type2 = 'T'; 
			else if(ordered_bps[i]->ty2 == C) type2 = 'C'; 
			else if(ordered_bps[i]->ty2 == G) type2 = 'G'; 
			//cout << type1 << " " << type2 << endl;
		}
			
		ordered_bps.clear();
		N = bps.size();
		counts++;
	}
	else {
		id = (int) compatible_bps.size()*udist01(rng);
		//cout << "Picked: " <<  indices[compatible_bps[id]] << endl;
		ordered_bps.push_back(bps[indices[compatible_bps[id]]]);
		for(int i = compatible_bps[id]; i < bps.size()-1; i++) {
			indices[i] = indices[i+1];
		}
		indices[bps.size()-1]=-1;
		for(int i = 0; i < bps.size(); i++) {
			//cout << " " << indices[i];
		}
		//cout << endl;
		N--;
	}	
}

for(int i = 0; i < ordered_bps.size(); i++) {
	if(ordered_bps[i]->ty1 == A) type1 = 'A';
	else if(ordered_bps[i]->ty1 == T) type1 = 'T'; 
	else if(ordered_bps[i]->ty1 == C) type1 = 'C'; 
	else if(ordered_bps[i]->ty1 == G) type1 = 'G';
	if(ordered_bps[i]->ty2 == A) type2 = 'A';
	else if(ordered_bps[i]->ty2 == T) type2 = 'T'; 
	else if(ordered_bps[i]->ty2 == C) type2 = 'C'; 
	else if(ordered_bps[i]->ty2 == G) type2 = 'G'; 
	cout << type1 << " " << type2 << endl;
}



//cout << "Generated!" << endl;

char type;
ofile.open("seq.txt");
ofile << "DOUBLE ";
for(int i = 0; i < ordered_bps.size(); i++) {
	if(ordered_bps[i]->ty1 == A) type = 'A';
	else if(ordered_bps[i]->ty1 == T) type = 'T'; 
	else if(ordered_bps[i]->ty1 == C) type = 'C'; 
	else if(ordered_bps[i]->ty1 == G) type = 'G'; 
	ofile << type;
}

if(ordered_bps[ordered_bps.size()-1]->ty2 == A) type = 'A';
else if(ordered_bps[ordered_bps.size()-1]->ty2 == T) type = 'T'; 
else if(ordered_bps[ordered_bps.size()-1]->ty2 == C) type = 'C'; 
else if(ordered_bps[ordered_bps.size()-1]->ty2 == G) type = 'G'; 
ofile << type;
ofile.close();

return 0;

}
