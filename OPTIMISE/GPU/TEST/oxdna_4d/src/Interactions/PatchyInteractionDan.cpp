/*
 * PatchyInteractionDan.cpp
 *
 *  Created on: 25/jan/2016
 *      Author: lorenzo -> dan
 */

//<iostream> is only included to allow cout, so I can print certain things more easily
#include <iostream>
#include "PatchyInteractionDan.h"
#include "../Utilities/Utils.h"

//Constructor

PatchyInteractionDan::PatchyInteractionDan() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(PATCHY, _patchy_interaction);

	/*//Set initialisation value to false
	 _initialised = false;*/

	//Set all arrays to NULL, so can delete them in destructor
	_N_particles_of_type = NULL;
	_particle_type_of = NULL;
	_N_patches_type = NULL;
	_patches_on_particle = NULL;
	_patch_vectors_type = NULL;
	_patch_vectors_particle = NULL;
	_patch_type_of = NULL;
	_sigma_ang_patch = NULL;
	_epsilon_patch = NULL;
	_ref_vectors_type = NULL;
	_ref_vectors_particle = NULL;
	_sigma_tor_patch = NULL;
	_number_offset_angles = NULL;
	_offset_angle_patch = NULL;
	/*particles = NULL;*/

}

//Destructor

PatchyInteractionDan::~PatchyInteractionDan() {
	//printf("PI, ~PatchyInteractionDan\n");

	//if(_initialised) {

	//Delete all arrays
	/*Template
	 if(_something != NULL) delete[] _something;
	 For pointers to pointers (etc.), you could do what is commented out for _patch_vectors_particle and _ref_vectors_particle, but that would be meaningful IFF all _patch_vectors_particle[i] were initialised to NULL just after you allocate memory for _patch_vectors_particle itself. i.e. In the code, after "_patch_vectors_particle = new LR_vector*[_N_particles];", I would have to loop through and set "_patch_vectors_particle[i] = NULL" for all i. There is no strong/immediate need to do this, so I haven't.*/
	if(_N_particles_of_type != NULL)
		delete[] _N_particles_of_type;
	if(_particle_type_of != NULL)
		delete[] _particle_type_of;
	if(_N_patches_type != NULL)
		delete[] _N_patches_type;
	if(_patches_on_particle != NULL)
		delete[] _patches_on_particle;
	if(_patch_vectors_type != NULL) {
		for(int part_type = 0; part_type < _N_particle_types; part_type++)
			delete[] _patch_vectors_type[part_type];
		delete[] _patch_vectors_type;
	}
	if(_patch_vectors_particle != NULL) {
		for(int particle = 0; particle < _N_particles; particle++)
			delete[] _patch_vectors_particle[particle];
		/*for (int particle = 0; particle < _N_particles; particle++) {
		 if(_patch_vectors_particle[particle] != NULL) delete[] _patch_vectors_particle[particle];
		 }*/
		delete[] _patch_vectors_particle;
	}
	if(_patch_type_of != NULL) {
		for(int part_type = 0; part_type < _N_particle_types; part_type++)
			delete[] _patch_type_of[part_type];
		delete[] _patch_type_of;
	}
	if(_sigma_ang_patch != NULL)
		delete[] _sigma_ang_patch;
	if(_epsilon_patch != NULL) {
		for(int patch_type = 0; patch_type < _N_patch_types; patch_type++)
			delete[] _epsilon_patch[patch_type];
		delete[] _epsilon_patch;
	}
	if(_ref_vectors_type != NULL) {
		for(int part_type = 0; part_type < _N_particle_types; part_type++)
			delete[] _ref_vectors_type[part_type];
		delete[] _ref_vectors_type;
	}
	if(_ref_vectors_particle != NULL) {
		for(int particle = 0; particle < _N_particles; particle++)
			delete[] _ref_vectors_particle[particle];
		/*for (int particle = 0; particle < _N_particles; particle++) {
		 if(_ref_vectors_particle[particle] != NULL) delete[] _ref_vectors_particle[particle];
		 }*/
		delete[] _ref_vectors_particle;
	}
	if(_sigma_tor_patch != NULL) {
		for(int patch_type = 0; patch_type < _N_patch_types; patch_type++)
			delete[] _sigma_tor_patch[patch_type];
		delete[] _sigma_tor_patch;
	}
	if(_number_offset_angles != NULL) {
		for(int patch_type = 0; patch_type < _N_patch_types; patch_type++)
			delete[] _number_offset_angles[patch_type];
		delete[] _number_offset_angles;
	}
	if(_offset_angle_patch != NULL) {
		for(int patch_type1 = 0; patch_type1 < _N_patch_types; patch_type1++) {
			for(int patch_type2 = 0; patch_type2 < _N_patch_types; patch_type2++)
				delete[] _offset_angle_patch[patch_type1][patch_type2];
			delete[] _offset_angle_patch[patch_type1];
		}
		delete[] _offset_angle_patch;
	}
	/*if(particles != NULL) {
	 for (int particle = 0; particle < _N_particles; particle++) delete[] particles[particle];
	 delete[] particles;
	 }*/

	//}
}

//Settings from input file

void PatchyInteractionDan::get_settings(input_file &inp) {
	//printf("PI, get_settings\n");

	//To temporarily store numbers (usually extracted from files)
	float tmp1;
	bool tmp2;

	//_rcut; code reused from original file
	tmp1 = DEFAULT_RCUT;
	getInputFloat(&inp, "PATCHY_rcut", &tmp1, 0);
	this->_rcut = (number) tmp1;

	//_tor_flag; code reused from original file
	tmp2 = DEFAULT_TOR_FLAG;
	getInputBool(&inp, "PATCHY_torsional", &tmp2, 0);
	_tor_flag = (bool) tmp2;

	//Gets topology filename
	BaseInteraction::get_settings(inp);

}

//Further initialisation

void PatchyInteractionDan::init() {
	//printf("PI, init\n");

	//Square of _rcut (for computational efficiency)
	this->_sqr_rcut = SQR(this->_rcut);

	//Repulsive LJ interaction at the cut-off
	//(epsilon would be * 4, sigma_LJ would be both numerators)
	_LJ_cut = 4 * (pow((1 / this->_sqr_rcut), 6) - pow((1 / this->_sqr_rcut), 3));

	//Calculate piecemeal potential crossover distance
	//(sigma_LJ would be in 2 numerators; 6 and 3 not 12 and 6 because use square of cutoff)
	/*number _r_crossover1 = (1 / this->_sqr_rcut);
	 printf("--> %.16lf\n", _r_crossover1);
	 number _r_crossover2 = pow((1 / this->_sqr_rcut), 6);
	 printf("--> %.16lf\n", _r_crossover2);
	 number _r_crossover3 = pow((1 / this->_sqr_rcut), 3);
	 printf("--> %.16lf\n", _r_crossover3);
	 number _r_crossover4 = -pow((1 / this->_sqr_rcut), 6) + pow((1 / this->_sqr_rcut), 3);
	 printf("--> %.16lf\n", _r_crossover4);
	 number _r_crossover5 = 4 * (-pow((1 / this->_sqr_rcut), 6) + pow((1 / this->_sqr_rcut), 3));
	 printf("--> %.16lf\n", _r_crossover5);
	 number _r_crossover6 = 1 - 4 * (-pow((1 / this->_sqr_rcut), 6) + pow((1 / this->_sqr_rcut), 3));
	 printf("--> %.16lf\n", _r_crossover6);
	 number _r_crossover7 = sqrt(1 - 4 * (-pow((1 / this->_sqr_rcut), 6) + pow((1 / this->_sqr_rcut), 3)) );
	 printf("--> %.16lf\n", _r_crossover7);
	 number _r_crossover8 = 1 + sqrt(1 - 4 * (-pow((1 / this->_sqr_rcut), 6) + pow((1 / this->_sqr_rcut), 3)) );
	 printf("--> %.16lf\n", _r_crossover8);
	 number _r_crossover9 = (1 + sqrt(1 - 4 * (-pow((1 / this->_sqr_rcut), 6) + pow((1 / this->_sqr_rcut), 3)) ) ) / 2;
	 printf("--> %.16lf\n", _r_crossover9);
	 number _r_crossover10 = pow(( (1 + sqrt(1 - 4 * (-pow((1 / this->_sqr_rcut), 6) + pow((1 / this->_sqr_rcut), 3)) ) ) / 2), (-1.0 / 6));
	 printf("--> %.16lf\n", _r_crossover10);*/
	_r_crossover = pow(((1 + sqrt(1 - 4 * (-pow((1 / this->_sqr_rcut), 6) + pow((1 / this->_sqr_rcut), 3)))) / 2), (-1.0 / 6));

	printf("IN1 TOLERANCE_PROJ %.16lf, _rcut %f, _sqr_rcut %f, _LJ_cut %f, _r_crossover %.16lf, _tor_flag %d\n", TOLERANCE_PROJ, this->_rcut, this->_sqr_rcut, _LJ_cut, _r_crossover, _tor_flag);
	//printf("IN1 TOLERANCE_PROJ %.16lf, TOLERANCE_PROD %.16lf, _rcut %f, _sqr_rcut %f, _LJ_cut %f, _r_crossover %.16lf, _tor_flag %d\n", TOLERANCE_PROJ, TOLERANCE_PROD, this->_rcut, this->_sqr_rcut, _LJ_cut, _r_crossover, _tor_flag);

	//count_pref = count_qref = count_cross = 0;

	OX_LOG(Logger::LOG_INFO, "Simulating a patchy (Dan) particle system.\n");

}

//Where is it called from, and so what are N, N_strands, etc.?

void PatchyInteractionDan::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
	//printf("PI, read_topology\n");

	//Reads lines from topology file into this
	char line[512];
	//Checks on numbers of particles, particle types, and patches 
	int particle_count_check = 0, part_type_count_check = 0, patch_count_check = 0;
	//To temporarily store numbers (usually extracted from files)
	double tmp1, tmp2, tmp3;

	//Tries to open topology file (filename obtained in get_settings)
	std::ifstream topology(this->_topology_filename, std::ios::in);
	if(!topology.good())
		throw oxDNAException("[In PatchyInteractionDan::get_settings] Can't read topology file '%s' . Aborting\n", this->_topology_filename);

	/*getInputString(&inp, "patch_file", this->_patch_filename, 1);
	 ifstream input_file;
	 input_file.open(patch_file);
	 if (!input_file.is_open())
	 {
	 printf("Patchy input file %s failed to open", patch_file);
	 filestr.close();
	 }*/

	//Until end of topology file
	while(!topology.eof()) {
		//printf("RT WHILE LOOP\n");

		//Read and parse first line
		topology.getline(line, 512);
		sscanf(line, "%d %d %d %*s\n", &_N_particles, &_N_particle_types, &_N_patch_types);
		printf("RT1 Number of particles %d, number of particle types %d, number of patch types %d\n", _N_particles, _N_particle_types, _N_patch_types);

		/*patch_file.getline(line, 512);
		 std::istringstream ss (line);
		 ss >> _N_particles >> _particle_types;*/

		//Create new arrays
		_N_particles_of_type = new int[_N_particle_types];
		_particle_type_of = new int[_N_particles];
		_N_patches_type = new int[_N_particle_types];
		_patches_on_particle = new int[_N_particles];
		//Define number of rows (one row for each particle type)
		_patch_vectors_type = new LR_vector*[_N_particle_types];
		//Define number of rows (one row for each particle)
		_patch_vectors_particle = new LR_vector*[_N_particles];
		//Define number of rows (one row for each particle type)
		_patch_type_of = new int*[_N_particle_types];
		_sigma_ang_patch = new number[_N_patch_types];
		//Define number of rows (one row for each patch type)
		_epsilon_patch = new number*[_N_patch_types];

		//Define number of rows (one row for each particle)
		_ref_vectors_particle = new LR_vector*[_N_particles];

		if(_tor_flag == true) {
			//Define number of rows (one row for each particle type)
			_ref_vectors_type = new LR_vector*[_N_particle_types];
			//Define number of rows (one row for each patch type)
			_sigma_tor_patch = new number*[_N_patch_types];
			//Define number of rows (one row for each patch type)
			_number_offset_angles = new int*[_N_patch_types];
			//Define number of rows (one row for each patch type)
			_offset_angle_patch = new number**[_N_patch_types];
		}

		//Iterate through particle types
		for(int part_type = 0; part_type < _N_particle_types; part_type++) {
			//printf("RT FOR LOOP\n");

			patch_count_check = 0;

			topology.getline(line, 512);
			/*istringstream ss(line);
			 ss >> _particles_of_type[part_type] >> _N_patches_type[part_type];*/
			sscanf(line, "%d %d %*s\n", &_N_particles_of_type[part_type], &_N_patches_type[part_type]);
			printf("RT2 part_type %d, number of particles of this type %d, number of patches on this type %d\n", part_type, _N_particles_of_type[part_type], _N_patches_type[part_type]);

			//Create new arrays
			//Define number of columns in each row (one row for each particle type, one column for each patch on that particle type)
			_patch_vectors_type[part_type] = new LR_vector[_N_patches_type[part_type]];
			//Define number of columns in each row (one row for each particle type, one column for each patch on that particle type)
			_patch_type_of[part_type] = new int[_N_patches_type[part_type]];

			//Iterate through patches on that particle type
			for(int patch = 0; patch < _N_patches_type[part_type]; patch++) {

				topology.getline(line, 512);
				sscanf(line, "%d %lf %lf %lf %*s\n", &_patch_type_of[part_type][patch], &tmp1, &tmp2, &tmp3);
				_patch_vectors_type[part_type][patch] = LR_vector(tmp1, tmp2, tmp3);
				printf("RT3a part_type %d, patch %d, patch_type %d, patchvect.x %.16lf, patchvect.y %.16lf patchvect.z %.16lf\n", part_type, patch, _patch_type_of[part_type][patch], _patch_vectors_type[part_type][patch].x, _patch_vectors_type[part_type][patch].y, _patch_vectors_type[part_type][patch].z);

				//printf("RT-- patchvect1 %.16lf, patchvect2 %.16lf patchvect3 %.16 refvect1 %.16lf, refvect2 %.16lf refvect3 %.16lf\n", tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);

				//printf(_patch_vectors_type[part_type][patch]);
				//std::cout << _patch_vectors_type[part_type][patch] << "\n";

				patch_count_check++;

			}

			//Check correct number of patches defined for this particle type [instance 1]
			if(patch_count_check != _N_patches_type[part_type])
				throw oxDNAException("Number of patches (%d) of particle type %d specified in topology file '%s' does not equal total number of patches (%d) specified above the list of patches [instance 1]. Aborting\n", patch_count_check, part_type, this->_topology_filename, _N_patches_type[part_type]);

			patch_count_check = 0;

			if(_tor_flag == true) {
				//Create new array
				//Define number of columns in each row (one row for each particle type, one column for each patch on that particle type)
				_ref_vectors_type[part_type] = new LR_vector[_N_patches_type[part_type]];

				//As above
				for(int patch = 0; patch < _N_patches_type[part_type]; patch++) {

					topology.getline(line, 512);
					sscanf(line, "%lf %lf %lf %*s\n", &tmp1, &tmp2, &tmp3);
					_ref_vectors_type[part_type][patch] = LR_vector(tmp1, tmp2, tmp3);
					printf("RT3b part_type %d, patch %d, patch_type %d, refvect.x %f, refvect.y %f refvect.z %f\n", part_type, patch, _patch_type_of[part_type][patch], _ref_vectors_type[part_type][patch].x, _ref_vectors_type[part_type][patch].y, _ref_vectors_type[part_type][patch].z);

					patch_count_check++;

				}

				//Check correct number of patches defined for this particle type [instance 2]
				if(patch_count_check != _N_patches_type[part_type])
					throw oxDNAException("Number of patches (%d) of particle type %d specified in topology file '%s' does not equal total number of patches (%d) specified above the list of patches [instance 2]. Aborting\n", patch_count_check, part_type, this->_topology_filename, _N_patches_type[part_type]);

			}

			part_type_count_check++;
			particle_count_check = particle_count_check + _N_particles_of_type[part_type];
			//printf("RT particle count check %d\n", particle_count_check);

		}

		/*topology.getline(line, 512);
		 sscanf(line, "%lf %*s\n", &_sigma_ang_patch[0]);//&tmp1
		 printf("*** _sigma_ang_patch[0] %f\n", _sigma_ang_patch[0]);//&tmp1*/

		//Check correct number of particles and particle types defined
		if(particle_count_check != _N_particles)
			throw oxDNAException("Sum of number of particles of each type (%d) in topology file '%s' does not equal total number of particles (%d) at top of file. Aborting\n", particle_count_check, this->_topology_filename, _N_particles);
		if(part_type_count_check != _N_particle_types)
			throw oxDNAException("Number of types of particles (%d) specified in topology file '%s' does not equal total number of particle types (%d) at top of file. Aborting\n", part_type_count_check, this->_topology_filename, _N_particle_types);

		//Iterate through patch types
		for(int patch_type = 0; patch_type < _N_patch_types; patch_type++) {

			topology >> _sigma_ang_patch[patch_type];
			topology.getline(line, 512);
			/*Gave warning message when compiling, because if use 'lf' it fails for float, and if use 'f' it fails for double.
			 topology.getline(line, 512);
			 sscanf(line, "%lf %*s\n", &_sigma_ang_patch[patch_type]);*/
			printf("RT4 _sigma_ang_patch[patch_type %d] %f\n", patch_type, _sigma_ang_patch[patch_type]);

			//Old code
			/*sscanf(line, "%lf %lf %lf %lf %*s\n", &_sigma_ang_patch[patch_type], &tmp1, &tmp2, &tmp3);
			 _ref_vector_patch[patch_type] = LR_vector(tmp1, tmp2, tmp3);
			 printf("RT4 patch_type %d, _sigma_ang_patch[patch_type] %f, _ref_vector_patch[patch_type] %lf %lf %lf\n", patch_type, _sigma_ang_patch[patch_type], _ref_vector_patch[patch_type].x, _ref_vector_patch[patch_type].y, _ref_vector_patch[patch_type].z);*/

		}

		/*FILE * top_file;
		 top_file = fopen(this->_topology_filename, "r");
		 //top_file.ignore(512,'\n');
		 fscanf(top_file, "%*s\n");

		 for (int patch_type = 0; patch_type < _N_patch_types; patch_type++) {
		 _sigma_tor_patch[patch_type] = new number[_N_patch_types];
		 //topology.getline(line, 512);
		 for (int patch_type2 = 0; patch_type2 < _N_patch_types; patch_type2++) {
		 fscanf(top_file, "%lf", &_sigma_tor_patch[patch_type][patch_type2]);
		 printf("RT5 _sigma_tor_patch[patch_type %d][patch_type2 %d] %f\n", patch_type, patch_type2, _sigma_tor_patch[patch_type][patch_type2]);
		 }
		 //sscanf(line, "%*s\n");
		 }

		 fclose(top_file);*/

		//Read in each row
		for(int patch_type1 = 0; patch_type1 < _N_patch_types; patch_type1++) {

			_epsilon_patch[patch_type1] = new number[_N_patch_types];

			std::cout << "RT5 _epsilon_patch[patch_type1 " << patch_type1 << "]";
			for(int patch_type2 = 0; patch_type2 < _N_patch_types; patch_type2++) {
				topology >> _epsilon_patch[patch_type1][patch_type2];
				std::cout << " " << _epsilon_patch[patch_type1][patch_type2];
			}
			std::cout << "\n";

			//Move to next row
			topology.getline(line, 512);

		}

		if(_tor_flag == true) {

			//Read in each row
			for(int patch_type1 = 0; patch_type1 < _N_patch_types; patch_type1++) {

				//Create new array
				//Define number of columns in each row (one row for each patch type, one column for each patch type)
				_sigma_tor_patch[patch_type1] = new number[_N_patch_types];

				std::cout << "RT6 _sigma_tor_patch[patch_type1 " << patch_type1 << "]";
				for(int patch_type2 = 0; patch_type2 < _N_patch_types; patch_type2++) {
					topology >> _sigma_tor_patch[patch_type1][patch_type2];
					std::cout << " " << _sigma_tor_patch[patch_type1][patch_type2];
				}
				std::cout << "\n";

				//Move to next row
				topology.getline(line, 512);

			}

			//Read in each row
			for(int patch_type1 = 0; patch_type1 < _N_patch_types; patch_type1++) {

				//Create new array
				//Define number of columns in each row (one row for each patch type, one column for each patch type)
				_number_offset_angles[patch_type1] = new int[_N_patch_types];
				_offset_angle_patch[patch_type1] = new number*[_N_patch_types];

				std::cout << "RT7 _offset_angle_patch[patch_type1 " << patch_type1 << "]";

				for(int patch_type2 = 0; patch_type2 < _N_patch_types; patch_type2++) {

					topology >> _number_offset_angles[patch_type1][patch_type2];

					//Create new array
					//Define number of values in each column in each row (one row for each patch type, one column for each patch type)
					_offset_angle_patch[patch_type1][patch_type2] = new number[_number_offset_angles[patch_type1][patch_type2]];

					std::cout << " " << _number_offset_angles[patch_type1][patch_type2] << ":";

					for(int offset_angle = 0; offset_angle < _number_offset_angles[patch_type1][patch_type2]; offset_angle++) {

						topology >> _offset_angle_patch[patch_type1][patch_type2][offset_angle];

						//[Revise this?] Check offset angle is in correct range
						if((_offset_angle_patch[patch_type1][patch_type2][offset_angle] > M_PI) || (_offset_angle_patch[patch_type1][patch_type2][offset_angle] <= -M_PI))
							throw oxDNAException("Offset angle %f, between patch type %d and patch type %d and offset angle number %d, specified in topology file '%s', is greater than PI or less than or equal to -PI. (Beware rounding of PI and -PI.) Aborting\n", _offset_angle_patch[patch_type1][patch_type2][offset_angle], patch_type1, patch_type2, offset_angle, this->_topology_filename);

						std::cout << " " << _offset_angle_patch[patch_type1][patch_type2][offset_angle] << ",";

					}

				}

				std::cout << "\n";

				//Move to next row
				topology.getline(line, 512);

			}

		}

		//Default values - redo?
		/*
		 //Sets default values of sigma_ang and sigma_tor; this block would usually not be within while loop
		 //DT - may reuse this block for sigma_LJ, epsilon, and other parameters
		 float tmp0 = 0.3;
		 //getInputFloat(&inp, "PATCHY_sigma", &tmp, 0);
		 _sigma_ang = (number) tmp0;
		 _sigma_tor = 2.0 * _sigma_ang;
		 //DT - old
		 //DTtmp = 0.12;
		 //DTgetInputFloat(&inp, "PATCHY_alpha", &tmp, 0);
		 //DT_patch_alpha = (number) tmp;*/

	}

	topology.close();

	/*OLD - early attempt
	 while (getline(patch_file, line)) {
	 istringstream ss(line);
	 patch_file >> _particle_types;
	 patch_file.ignore(256, '\n');
	 for (int part_type = 0; part_type < _particle_types; part_type++) {
	 patch_file >> _particles_of_type[part_type] >> _N_patches_type[part_type];
	 patch_file.ignore(256, '\n');
	 }
	 }
	 this->N_int_centers() = N_patches;

	 input_file.close();*/

	//N_strands is used elsewhere in the code
	*N_strands = _N_particles;

	//This needs to be called here (because it uses 'particles')
	allocate_particles(particles);

	//Define particle index and type
	for(int particle_number = 0; particle_number < _N_particles; particle_number++) {
		particles[particle_number]->index = particle_number;
		//Could move this to PatchyParticle.cpp in the future
		particles[particle_number]->type = _particle_type_of[particle_number];
		//strand_id is used elsewhere in the code
		particles[particle_number]->strand_id = particle_number;
	}

	/*printf("COMPARISON WITH EVA\n");
	 printf("p->index, q->index, r->x, r->y, r->z, sqr_r_dist, r_dist, V_LJ (before shift), V_LJ (after shift)\n");
	 printf(", , p_patch, q_patch, ppatch.x, ppatch.y, ppatch.z, qpatch.x, qpatch.y, qpatch.z, pref.x, pref.y, pref.z, qref.x, qref.y, qref.z, angle_r_p, angle_r_q, V_ang_p, V_ang_q, V_ang, proj_pref.x, proj_pref.y, proj_pref.z, proj_qref.x, proj_qref.y, proj_qref.z, angle_tor (initial), cross_proj </>, dot_proj </>, sign </>, angle_tor (final), offset_angle1, angle_diff1, offset_angle2, angle_diff2, offset_angle3, angle_diff3, min_sqr_angle_diff, V_tor, V_ang_V_tor_epsilon, max_V_ang_V_tor_epsilon (at current stage)\n");
	 printf(", max_V_ang_V_tor_epsilon (final), energy\n");*/

	/*//Set initialisation flag to true, now all pointers are initialised
	 _initialised = true;*/

}

/*OLD 3/6/16

 void PatchyInteractionDan::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
 //printf("PI, read_topology\n");

 //Is this needed?
 *N_strands = N;

 //Does this need to be called here?
 allocate_particles(particles, N);

 //Can I safely delete all this?
 for (int i = 0; i < N; i ++) {
 //particles[i]->type = (i < _N_A) ? P_A : P_B;
 //particles[i]->btype = (i < _N_A) ? P_A : P_B;
 particles[i]->strand_id = i;
 }

 }*/

void PatchyInteractionDan::allocate_particles(std::vector<BaseParticle*> &particles) {
	//printf("PI, allocate_particles\n");

	int particle_number = 0;

	//Define properties of individual particles, based on information read in read_topology
	for(int part_type = 0; part_type < _N_particle_types; part_type++) {

		for(int type_particle_count = 0; type_particle_count < _N_particles_of_type[part_type]; type_particle_count++) {

			//Particle type and number of patches
			_particle_type_of[particle_number] = part_type;
			_patches_on_particle[particle_number] = _N_patches_type[part_type];

			//Create new arrays
			//Define number of columns in each row (one row for each particle, one column for each patch on that particle)
			_patch_vectors_particle[particle_number] = new LR_vector[_N_patches_type[part_type]];
			_ref_vectors_particle[particle_number] = new LR_vector[_N_patches_type[part_type]];

			for(int patch = 0; patch < _N_patches_type[part_type]; patch++) {
				//printf("AP0 _patch_type_of[part_type %d][patch %d] %d\n", part_type, patch, _patch_type_of[part_type][patch]);

				_patch_vectors_particle[particle_number][patch] = _patch_vectors_type[part_type][patch];
				if(_tor_flag == true) {
					_ref_vectors_particle[particle_number][patch] = _ref_vectors_type[part_type][patch];
				}
				else {
					_ref_vectors_particle[particle_number][patch] = _patch_vectors_particle[particle_number][patch];
				}

				//printf("AP1 _N_patches_type[part_type] %d, _patch_vectors_particle[particle_number][patch] %f %f %f, _ref_vectors_particle[particle_number][patch] %f %f %f\n", _N_patches_type[part_type], _patch_vectors_particle[particle_number][patch].x, _patch_vectors_particle[particle_number][patch].y, _patch_vectors_particle[particle_number][patch].z, _ref_vectors_particle[particle_number][patch].x, _ref_vectors_particle[particle_number][patch].y, _ref_vectors_particle[particle_number][patch].z);

			}

			//printf("AP2 part_type %d, type_particle_count %d, particle_number %d, _particle_type_of[particle_number] %d, _patches_on_particle[particle_number] %d, _N_particles_of_type[part_type] %d\n", part_type, type_particle_count, particle_number, _particle_type_of[particle_number], _patches_on_particle[particle_number], _N_particles_of_type[part_type]);

			particle_number = particle_number + 1;

		}

	}

	//Check correct number of particles defined
	if(particle_number != _N_particles)
		throw oxDNAException("Number of particles allocated (%d) in 'allocate_particles' does not equal total number of particles %d. Aborting\n", particle_number, _N_particles);

	/*int part_type = 0;
	 int type_particle_count = 0;

	 for (int particle_number = 0; particle_number < _N_particles; particle_number++) {
	 _particle_type_of[particle_number] = part_type;
	 _patches_on_particle[particle_number] = _N_patches_type[part_type];
	 printf("AP part_type %d, type_particle_count %d, particle_number %d, _particle_type_of[particle_number] %d, _patches_on_particle[particle_number] %d, _N_particles_of_type[part_type] %d\n", part_type, type_particle_count, particle_number, _particle_type_of[particle_number], _patches_on_particle[particle_number], _N_particles_of_type[part_type]);
	 type_particle_count = type_particle_count + 1;
	 if (type_particle_count == _N_particles_of_type[part_type]) {
	 printf("type_particle_count %d MATCHES _N_particles_of_type[part_type] %d\n", type_particle_count, _N_particles_of_type[part_type]);
	 part_type = part_type + 1;
	 type_particle_count = 0;
	 }
	 }*/

	//Create particles
	for(particle_number = 0; particle_number < _N_particles; particle_number++) {
		//printf("Allocating particle %d, _patches_on_particle %d, _particle_type_of %d\n", particle_number, _patches_on_particle[particle_number], _particle_type_of[particle_number]);

		particles[particle_number] = new PatchyParticleDan(_patches_on_particle[particle_number], _patch_vectors_particle[particle_number], _ref_vectors_particle[particle_number], _tor_flag);

		//printf("Allocated particle %d\n", particle_number);

	}

}

//All interactions are nonbonded
number PatchyInteractionDan::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	//printf("PI, pair_interaction\n");

	//Currently only MC (not MD) simulations supported
	if(update_forces) {
		throw oxDNAException("PatchyInteractionDan does not support the calculation of forces and torques. Aborting\n");
	}

	return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

//No bonded interaction, so always 0
number PatchyInteractionDan::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return (number) 0.f;
}

number PatchyInteractionDan::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	LR_vector computed_r(0, 0, 0);

	//?? If r is not given, it computes it itself
	if(compute_r) {
		_computed_r = this->_box->min_image(p->pos, q->pos);
	}

	return _patchy_interaction(p, q, false, update_forces);
}

void PatchyInteractionDan::check_input_sanity(std::vector<BaseParticle*> &particles) {

}
