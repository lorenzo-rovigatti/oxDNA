/*
 * utils.h
 *
 *  Created on: Feb 19, 2013
 *      Author: petr
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

void load_fitting_options(const char *initfile, vector<string>& sequences, vector<string>& datafiles)
{
 ifstream fin(initfile);
 if (!fin.good())
 {
	 throw string(string("Error opening file ") + string(initfile));
 }

 string sequence_filename;
 string parameter_filename;

 fin >> sequence_filename >> parameter_filename;
 while(fin.good())
 {
	 sequences.push_back(sequence_filename);
	 datafiles.push_back(parameter_filename);
	 fin >> sequence_filename >> parameter_filename;
 }
}


#endif /* UTILS_H_ */
