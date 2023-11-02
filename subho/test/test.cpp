#include <iostream>
#include "Analysis.h"
using namespace std;

int main(){
    Analysis ico("ico.top","ico.dat","DNA"); // Loading the top and the dat files
    // ico.shiftbox({0,0,0}); // Bringing everything to 1st coordinate
    // ico.testBoxOverloaded(); // Checking overload from previous action.

    vector<vector<int>> ids{
        {6532,6555,6559,6589,6158,6217,4906,6182,6180,4925},
        {4471,4492,4497,5028,5033,5054,5059,5080,5085,4466},
        {4023,4044,4050,4070,6013,6033,6359,6380,6406,6386},
        {2205,6072,6077,6099,5232,5211,2515,6512,5257,2826},
        {2566,2781,3668,3497,3522,3549,3529,3554,3637,2722},
        {4678,1378,4755,4705,4781,4762,1484,1529,4672,4814},
        {4352,4179,4186,4207,4213,4294,4300,4320,4326,2004},
        {5389,746,695,128,5422,5442,5306,5449,5312,5333},
        {3186,3207,3233,3213,3239,2885,2877,3378,3384,3180},
        {6600,6620,275,6646,6314,6293,390,436,5735,5714},
        {2360,2318,2311,2264,2256,2471,2421,2461,3773,3794},
        {6225,6246,224,180,171,652,6665,6687,589,6712}
        };
    
    vector<int> colors{20,-20,30,-30,20,-20,35,-35,40,-40,35,-35}; // if 12 colors provided 13th color will be by default 100 or colorless
    vector<double> radius{8,10}; // if 2 of the radius are provided then 1st for all particle and 2nd for the central particle
    // cout << ids.size()<<endl;
    Analysis psp("","","newPSP"); // generating empty particle class
    ico.generatePSP(&psp,ids,colors,radius,5,10); // generating ccg particle from oxDNA file.
    psp.populate(125,1); //generate 125 crystal with extra 1 su seperation between them
    psp.boxToCubic();// make the box cubic nothing is effected
    psp.writeCCGtopology("125.top"); // write the CCG topology
    psp.writeCCGviewTopology("125iew.top"); // write oxview readeable topology
    psp.writeConfig("125.dat"); // write config dat file
    // cout <<psp.particles[5].connector<<endl;
    return 0;
}