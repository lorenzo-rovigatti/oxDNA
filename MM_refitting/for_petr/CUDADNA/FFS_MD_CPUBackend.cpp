/*
 * MD_CPUBackend.cpp
 *
 *  Created on: 03/set/2010
 *      Author: sulc
 */

//#include "SimManager.h"
#include "FFS_MD_CPUBackend.h"
#include "IOManager.h"
#include <sstream>
#include <cctype>

#include "SimManager.h"
//extern bool SimManager::stop;

using namespace std;

bool parsed_expression::parse_expression(const char *expression,
		OrderParameters *op) {

	string s(expression);
	//s.erase(remove_if(s.begin(), s.end(),' '), s.end()); //remove spaces
	s.erase(remove_if(s.begin(), s.end(), static_cast<int(*)(int)>( isspace )), s.end());

	//cout << "Processing " << s << endl;
	int condition_type = -1;
	int pos = 0;
	unsigned int i;
	for (i = 0; i < s.length() - 1; i++) {
		if (s[i] == '<') {
			if (s[i + 1] == '=') {
				condition_type = LEQ;
				pos = i + 2;
				break;
			} else {
				condition_type = LT;
				pos = i + 1;
				break;
			}
		} else if (s[i] == '>') {
			if (s[i + 1] == '=') {
				condition_type = GEQ;
				pos = i + 2;
				break;
			} else {
				condition_type = GT;
				pos = i + 1;
				break;
			}

		} else if (s[i] == '=' && s[i + 1] == '=') {
			condition_type = EQ;
			pos = i + 2;
			break;
		}

	}
	if (condition_type == -1) {
		return false;
	}
	compare_type = condition_type;

	//cout << "Processed as " << condition_type << " and left with  "
	//		<< s.substr(0, i) << endl;

	string parname = s.substr(0, i);
	parameter_index = op->get_hbpar_id_from_name(parname.c_str());
	expression_type = 0;
	if (parameter_index == -1) { //it is a distance parameter

		parameter_index = op->get_distpar_id_from_name(parname.c_str());
		expression_type = 1;
	}
	if (parameter_index == -1) {
		//cout << "did not find parameter " << parname.c_str() << endl;
		return false;
	}

	stringstream str;
	str << s.substr(pos);

	if (expression_type == 0)
		str >> value_to_compare_with_hb;
	else if (expression_type == 1)
		str >> value_to_compare_with_dist;

	return true;

}

bool parsed_expression::eval_expression(OrderParameters *op) {

	if (expression_type == 0) {

		int current_value = op->get_hb_parameter(parameter_index);
		//cout << "Evaluating HB with " << current_value << " to be " << expression_type << " with " << value_to_compare_with_hb << endl;
		switch (compare_type) {
		case GEQ:
			return current_value >= value_to_compare_with_hb ? true : false;
			break;

		case LEQ:
			return current_value <= value_to_compare_with_hb ? true : false;
			break;

		case LT:
			return current_value < value_to_compare_with_hb ? true : false;
			break;
		case GT:
			return current_value > value_to_compare_with_hb ? true : false;
			break;
		case EQ:
			return current_value == value_to_compare_with_hb ? true : false;
			break;
		default:
			break;
		}
	} else if (expression_type == 1) {
		double current_value = op->get_distance_parameter(parameter_index);
		switch (expression_type) {
		case GEQ:
			return current_value >= value_to_compare_with_dist ? true : false;
			break;

		case LEQ:
			return current_value <= value_to_compare_with_dist ? true : false;
			break;

		case LT:
			return current_value < value_to_compare_with_dist ? true : false;
			break;
		case GT:
			return current_value > value_to_compare_with_dist ? true : false;
			break;
		case EQ:
			return current_value == value_to_compare_with_dist ? true : false;
			break;
		default:
			break;
		}
	} else {
		return false;
	}

	return false;

}

template<typename number>
FFS_MD_CPUBackend<number>::FFS_MD_CPUBackend(IOManager *IO) :
		MD_CPUBackend<number>(IO) {
	this->_is_CUDA_sim = false;
	_ffs_type = -1;
}

template<typename number>
FFS_MD_CPUBackend<number>::~FFS_MD_CPUBackend() {

}

template<typename number>
inline number FFS_MD_CPUBackend<number>::_ffs_particle_particle_interaction(
		Particle<number> *p, Particle<number> *q) {
	// true if p and q are Watson-Crick pairs
	//bool is_pair = (q->type + p->type == 3);
	bool is_pair = (q->btype + p->btype == 3);

	LR_vector<number> r = q->pos.minimum_image(p->pos, this->_box_side);

	number energy = 0;

	// excluded volume

	// BASE-BASE
	LR_vector<number> force(0, 0, 0);
	LR_vector<number> rcenter = r + q->pos_base - p->pos_base;
	energy += _excluded_volume(rcenter, force, EXCL_S2, EXCL_R2, EXCL_B2,
			EXCL_RC2);
	LR_vector<number> torquep = -p->pos_base.cross(force);
	LR_vector<number> torqueq = q->pos_base.cross(force);

	p->force -= force;
	q->force += force;

	// P-BASE vs. Q-BACK
	rcenter = r + q->pos_back - p->pos_base;
	energy += _excluded_volume(rcenter, force, EXCL_S3, EXCL_R3, EXCL_B3,
			EXCL_RC3);
	torquep += -p->pos_base.cross(force);
	torqueq += q->pos_back.cross(force);

	p->force -= force;
	q->force += force;

	// P-BACK vs. Q-BASE
	rcenter = r + q->pos_base - p->pos_back;
	energy += _excluded_volume(rcenter, force, EXCL_S4, EXCL_R4, EXCL_B4,
			EXCL_RC4);
	torquep += -p->pos_back.cross(force);
	torqueq += q->pos_base.cross(force);

	p->force -= force;
	q->force += force;

	// BACK-BACK
	rcenter = r + q->pos_back - p->pos_back;
	energy += _excluded_volume(rcenter, force, EXCL_S1, EXCL_R1, EXCL_B1,
			EXCL_RC1);
	torquep += -p->pos_back.cross(force);
	torqueq += q->pos_back.cross(force);

	p->force -= force;
	q->force += force;

	// HYDROGEN BONDING
	LR_vector<number> rhydro = r + q->pos_base - p->pos_base;
	number rhydromod = rhydro.module();
	if (is_pair && HYDR_RCLOW < rhydromod && rhydromod < HYDR_RCHIGH) {
		// vector, versor and magnitude of the base-base separation
		LR_vector<number> rhydrodir = rhydro / rhydromod;

		// particle axes according to Allen's paper
		LR_vector<number> a1 = p->orientationT.v1;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> b3 = q->orientationT.v3;

		// angles involved in the HB interaction
		/*number t1 = LRACOS (-a1 * b1);
		 number t2 = LRACOS (-b1 * rhydrodir);
		 number t3 = LRACOS ( a1 * rhydrodir);
		 number t4 = LRACOS ( a3 * b3);
		 number t7 = LRACOS (-b3 * rhydrodir);
		 number t8 = LRACOS ( a3 * rhydrodir);
		 */
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rhydrodir;
		number cost3 = a1 * rhydrodir;

		number cost4 = a3 * b3;
		number cost7 = -b3 * rhydrodir;
		number cost8 = a3 * rhydrodir;

		// functions called at their relevant arguments
		number f1 = this->_interaction.f1(rhydromod, HYDR_F1, q->type, p->type);
		/*
		 number f4t1 = this->_interaction.f4(t1, HYDR_F4_THETA1);
		 number f4t2 = this->_interaction.f4(t2, HYDR_F4_THETA2);
		 number f4t3 = this->_interaction.f4(t3, HYDR_F4_THETA3);
		 */
		number f4t1 = this->_interaction.query_mesh(cost1,
				this->_interaction.mesh_f4[HYDR_F4_THETA1]);
		number f4t2 = this->_interaction.query_mesh(cost2,
				this->_interaction.mesh_f4[HYDR_F4_THETA2]);
		number f4t3 = this->_interaction.query_mesh(cost3,
				this->_interaction.mesh_f4[HYDR_F4_THETA3]);

		number f4t4 = this->_interaction.query_mesh(cost4,
				this->_interaction.mesh_f4[HYDR_F4_THETA4]);
		number f4t7 = this->_interaction.query_mesh(cost7,
				this->_interaction.mesh_f4[HYDR_F4_THETA7]);
		number f4t8 = this->_interaction.query_mesh(cost8,
				this->_interaction.mesh_f4[HYDR_F4_THETA8]);
		/*number f4t4 = this->_interaction.f4(t4, HYDR_F4_THETA4);
		 number f4t7 = this->_interaction.f4(t7, HYDR_F4_THETA7);
		 number f4t8 = this->_interaction.f4(t8, HYDR_F4_THETA8);*/

		number hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		if (hb_energy < BOND_CUTOFF) {
			_op.add_hb(q->index, p->index);

		}

		energy += hb_energy;

		this->_U_hydr += hb_energy;

		// makes sense, since the above functions may return 0. exactly
		if (hb_energy != 0.) {
			// derivatives called at the relevant arguments
			number f1D = this->_interaction.f1D(rhydromod, HYDR_F1, q->type,
					p->type);
			/*number f4t1Dsin = -this->_interaction.f4Dsin(t1, HYDR_F4_THETA1);
			 number f4t2Dsin = -this->_interaction.f4Dsin(t2, HYDR_F4_THETA2);
			 number f4t3Dsin =  this->_interaction.f4Dsin(t3, HYDR_F4_THETA3);*/
			number f4t1Dsin = this->_interaction.query_meshD(cost1,
					this->_interaction.mesh_f4[HYDR_F4_THETA1]);
			number f4t2Dsin = this->_interaction.query_meshD(cost2,
					this->_interaction.mesh_f4[HYDR_F4_THETA2]);
			number f4t3Dsin = -this->_interaction.query_meshD(cost3,
					this->_interaction.mesh_f4[HYDR_F4_THETA3]);

			number f4t4Dsin = -this->_interaction.query_meshD(cost4,
					this->_interaction.mesh_f4[HYDR_F4_THETA4]);
			number f4t7Dsin = this->_interaction.query_meshD(cost7,
					this->_interaction.mesh_f4[HYDR_F4_THETA7]);
			number f4t8Dsin = -this->_interaction.query_meshD(cost8,
					this->_interaction.mesh_f4[HYDR_F4_THETA8]);
			/*number f4t4Dsin =  this->_interaction.f4Dsin(t4, HYDR_F4_THETA4);
			 number f4t7Dsin = -this->_interaction.f4Dsin(t7, HYDR_F4_THETA7);
			 number f4t8Dsin =  this->_interaction.f4Dsin(t8, HYDR_F4_THETA8);*/

			//f4t1 = 1.;
			//f4t2 = 1.;
			//f4t3 = 1.;
			//f4t4 = 1.;
			//f4t7 = 1.;
			//f4t8 = 1.;
			//f4t1Dsin = 0.;
			//f4t2Dsin = 0.;
			//f4t3Dsin = 0.;
			//f4t4Dsin = 0.;
			//f4t7Dsin = 0.;
			//f4t8Dsin = 0.;
			// RADIAL PART
			force = -rhydrodir * f1D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

			// TETA4; t4 = LRACOS (a3 * b3);
			LR_vector<number> dir = a3.cross(b3);
			number torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7
					* f4t8;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// TETA1; t1 = LRACOS (-a1 * b1);
			dir = a1.cross(b1);
			torquemod = -f1 * f4t1Dsin * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// TETA2; t2 = LRACOS (-b1 * rhydrodir);
			number fact = f1 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;
			//force += (b1 + rhydrodir * cos(t2)) / rhydromod * fact;
			force += (b1 + rhydrodir * cost2) / rhydromod * fact;
			dir = rhydrodir.cross(b1);
			//torquemod = - f1 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;

			torqueq += -dir * fact;

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
			//force += (a1 - rhydrodir * cos(t3)) / rhydromod * f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;
			force += (a1 - rhydrodir * cost3) / rhydromod * f1 * f4t1 * f4t2
					* f4t3Dsin * f4t4 * f4t7 * f4t8;

			LR_vector<number> t3dir = rhydrodir.cross(a1);
			torquemod = -f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			torquep += t3dir * torquemod;

			// THETA7; t7 = LRACOS (-rhydrodir * b3);
			//force += (b3 + rhydrodir * cos(t7)) / rhydromod * f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;
			force += (b3 + rhydrodir * cost7) / rhydromod * f1 * f4t1 * f4t2
					* f4t3 * f4t4 * f4t7Dsin * f4t8;

			LR_vector<number> t7dir = rhydrodir.cross(b3);
			torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			torqueq += t7dir * torquemod;

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			//force +=  (a3 - rhydrodir * cos(t8)) / rhydromod * f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;
			force += (a3 - rhydrodir * cost8) / rhydromod * f1 * f4t1 * f4t2
					* f4t3 * f4t4 * f4t7 * f4t8Dsin;

			LR_vector<number> t8dir = rhydrodir.cross(a3);
			torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			torquep += t8dir * torquemod;

			// final update of forces and torques for HB
			p->force -= force;
			q->force += force;

			torquep -= p->pos_base.cross(force);
			torqueq += q->pos_base.cross(force);
		}
	}
	// END OF HYDROGEN BONDING

	// CROSS STACKING
	LR_vector<number> rcstack = rhydro;
	//LR_vector<number> rcstack = r + q->pos_base - p->pos_base;
	number rcstackmod = rhydromod;
	//number rcstackmod = rcstack.module();
	if (CRST_RCLOW < rcstackmod && rcstackmod < CRST_RCHIGH) {
		LR_vector<number> rcstackdir = rcstack / rcstackmod;

		// particle axes according to Allen's paper
		LR_vector<number> a1 = p->orientationT.v1;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> b3 = q->orientationT.v3;

		// angles involved in the CRST interaction
		/*
		 number t1 = LRACOS (-a1 * b1);
		 number t2 = LRACOS (-b1 * rcstackdir);
		 number t4 = LRACOS ( a3 * b3);
		 number t3 = LRACOS ( a1 * rcstackdir);
		 number t7 = LRACOS (-rcstackdir * b3);
		 number t8 = LRACOS ( rcstackdir * a3);
		 */
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rcstackdir;
		number cost3 = a1 * rcstackdir;
		number cost4 = a3 * b3;
		number cost7 = -b3 * rcstackdir;
		number cost8 = a3 * rcstackdir;

		// functions called at their relevant arguments
		number f2 = this->_interaction.f2(rcstackmod, CRST_F2);
		/*
		 number f4t1 = this->_interaction.f4(t1, CRST_F4_THETA1);
		 number f4t2 = this->_interaction.f4(t2, CRST_F4_THETA2);
		 number f4t3 = this->_interaction.f4(t3, CRST_F4_THETA3);
		 number f4t4 = this->_interaction.f4(t4, CRST_F4_THETA4) + this->_interaction.f4(PI - t4, CRST_F4_THETA4);
		 */
		number f4t1 = this->_interaction.query_mesh(cost1,
				this->_interaction.mesh_f4[CRST_F4_THETA1]);
		number f4t2 = this->_interaction.query_mesh(cost2,
				this->_interaction.mesh_f4[CRST_F4_THETA2]);
		number f4t3 = this->_interaction.query_mesh(cost3,
				this->_interaction.mesh_f4[CRST_F4_THETA3]);
		number f4t4 = this->_interaction.query_mesh(cost4,
				this->_interaction.mesh_f4[CRST_F4_THETA4])
				+ this->_interaction.query_mesh(-cost4,
						this->_interaction.mesh_f4[CRST_F4_THETA4]);
		/*
		 number f4t7 = this->_interaction.f4(t7, CRST_F4_THETA7) + this->_interaction.f4(PI - t7, CRST_F4_THETA7);
		 number f4t8 = this->_interaction.f4(t8, CRST_F4_THETA8) + this->_interaction.f4(PI - t8, CRST_F4_THETA8);
		 */
		number f4t7 = this->_interaction.query_mesh(cost7,
				this->_interaction.mesh_f4[CRST_F4_THETA7])
				+ this->_interaction.query_mesh(-cost7,
						this->_interaction.mesh_f4[CRST_F4_THETA7]);
		;
		number f4t8 = this->_interaction.query_mesh(cost8,
				this->_interaction.mesh_f4[CRST_F4_THETA8])
				+ this->_interaction.query_mesh(-cost8,
						this->_interaction.mesh_f4[CRST_F4_THETA8]);
		;

		number cstk_energy = f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
		energy += cstk_energy;

		// makes sense since the above f? can return exacly 0.
		if (cstk_energy != 0.) {
			// derivatives called at the relevant arguments
			number f2D = this->_interaction.f2D(rcstackmod, CRST_F2);
			/*number f4t1Dsin = -this->_interaction.f4Dsin(t1, CRST_F4_THETA1);
			 number f4t2Dsin = -this->_interaction.f4Dsin(t2, CRST_F4_THETA2);
			 number f4t3Dsin =  this->_interaction.f4Dsin(t3, CRST_F4_THETA3);
			 number f4t4Dsin =  this->_interaction.f4Dsin(t4, CRST_F4_THETA4) - this->_interaction.f4Dsin(PI - t4, CRST_F4_THETA4);*/
			number f4t1Dsin = this->_interaction.query_meshD(cost1,
					this->_interaction.mesh_f4[CRST_F4_THETA1]);
			number f4t2Dsin = this->_interaction.query_meshD(cost2,
					this->_interaction.mesh_f4[CRST_F4_THETA2]);
			number f4t3Dsin = -this->_interaction.query_meshD(cost3,
					this->_interaction.mesh_f4[CRST_F4_THETA3]);
			number f4t4Dsin = -this->_interaction.query_meshD(cost4,
					this->_interaction.mesh_f4[CRST_F4_THETA4])
					+ this->_interaction.query_meshD(-cost4,
							this->_interaction.mesh_f4[CRST_F4_THETA4]);
			/*
			 number f4t7Dsin = -this->_interaction.f4Dsin(t7, CRST_F4_THETA7) + this->_interaction.f4Dsin(PI - t7, CRST_F4_THETA7);
			 number f4t8Dsin =  this->_interaction.f4Dsin(t8, CRST_F4_THETA8) - this->_interaction.f4Dsin(PI - t8, CRST_F4_THETA8);
			 */
			number f4t7Dsin = this->_interaction.query_meshD(cost7,
					this->_interaction.mesh_f4[CRST_F4_THETA7])
					- this->_interaction.query_meshD(-cost7,
							this->_interaction.mesh_f4[CRST_F4_THETA7]);
			number f4t8Dsin = -this->_interaction.query_meshD(cost8,
					this->_interaction.mesh_f4[CRST_F4_THETA8])
					+ this->_interaction.query_meshD(-cost8,
							this->_interaction.mesh_f4[CRST_F4_THETA8]);

			//f4t1 = 1.;
			//f4t2 = 1.;
			//f4t3 = 1.;
			//f4t4 = 1.;
			//f4t7 = 1.;
			//f4t8 = 1.;

			//f4t1Dsin = 0.;
			//f4t2Dsin = 0.;
			//f4t3Dsin = 0.;
			//f4t4Dsin = 0.;
			//f4t7Dsin = 0.;
			//f4t8Dsin = 0.;

			/*
			 FILE * yo = fopen("yo.dat", "w");
			 for (int kk = 0; kk<1000.; kk++)
			 {
			 fprintf(yo,"%lf %lf %lf\n", 0.002 * (kk - 500) * 2 * PI, this->_interaction.f4(0.002 * (kk-500) * 2 * PI, CRST_F4_THETA3), this->_interaction.f4D(0.002 * (kk-500) * 2 * PI, CRST_F4_THETA3));
			 }
			 fclose(yo);
			 exit(-1);*/

			// RADIAL PART
			force = -rcstackdir * f2D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

			// THETA1; t1 = LRACOS (-a1 * b1);
			LR_vector<number> t1dir = a1.cross(b1);
			number torquemod = -f2 * f4t1Dsin * f4t2 * f4t3 * f4t4 * f4t7
					* f4t8;

			torquep -= t1dir * torquemod;
			torqueq += t1dir * torquemod;

			// TETA2; t2 = LRACOS (-b1 * rhydrodir);
			//force += (b1 + rcstackdir * cos(t2)) / rcstackmod * f2 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;
			force += (b1 + rcstackdir * cost2) / rcstackmod * f2 * f4t1
					* f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;

			LR_vector<number> t2dir = rcstackdir.cross(b1);
			torquemod = -f2 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;

			torqueq += t2dir * torquemod;

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
			//force += (a1 - rcstackdir * cos(t3)) / rcstackmod * f2 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;
			force += (a1 - rcstackdir * cost3) / rcstackmod * f2 * f4t1 * f4t2
					* f4t3Dsin * f4t4 * f4t7 * f4t8;

			LR_vector<number> t3dir = rcstackdir.cross(a1);
			torquemod = -f2 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			torquep += t3dir * torquemod;

			// TETA4; t4 = LRACOS (a3 * b3);
			LR_vector<number> t4dir = a3.cross(b3);
			torquemod = -f2 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7 * f4t8;

			torquep -= t4dir * torquemod;
			torqueq += t4dir * torquemod;

			// THETA7; t7 = LRACOS (-rcsrackir * b3);
			//force += (b3 + rcstackdir * cos(t7)) / rcstackmod * f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;
			force += (b3 + rcstackdir * cost7) / rcstackmod * f2 * f4t1 * f4t2
					* f4t3 * f4t4 * f4t7Dsin * f4t8;

			LR_vector<number> t7dir = rcstackdir.cross(b3);
			torquemod = -f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			torqueq += t7dir * torquemod;

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			// force +=  (a3 - rcstackdir * cos(t8)) / rcstackmod * f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;
			force += (a3 - rcstackdir * cost8) / rcstackmod * f2 * f4t1 * f4t2
					* f4t3 * f4t4 * f4t7 * f4t8Dsin;

			LR_vector<number> t8dir = rcstackdir.cross(a3);
			torquemod = -f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			torquep += t8dir * torquemod;

			// final update of forces and torques for CRST
			p->force -= force;
			q->force += force;

			torquep -= p->pos_base.cross(force);
			torqueq += q->pos_base.cross(force);
		}
	}

	// COAXIAL STACKING
	LR_vector<number> rstack = r + q->pos_stack - p->pos_stack;
	number rstackmod = rstack.module();
	if (CXST_RCLOW < rstackmod && rstackmod < CXST_RCHIGH) {
		LR_vector<number> rstackdir = rstack / rstackmod;

		// particle axes according to Allen's paper
		LR_vector<number> a1 = p->orientationT.v1;
		LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> b3 = q->orientationT.v3;

		// angles involved in the CXST interaction
		/*
		 number t1 = LRACOS (-a1 * b1);
		 number t4 = LRACOS ( a3 * b3);
		 number t5 = LRACOS ( a3 * rstackdir);
		 number t6 = LRACOS (-b3 * rstackdir);
		 */
		number cost1 = -a1 * b1;
		number cost4 = a3 * b3;
		number cost5 = a3 * rstackdir;
		number cost6 = -b3 * rstackdir;
		// here we define rbackghost, which is the position the backbone would be in if we were using the no-grooving model. We do this because some of the force calculations rely on this definition of rback; it's just easier than recalculating the constants in that calculation
		// NB if major-minor grooving is not in use, rback = rbackboneghost and everything works as it should
		LR_vector<number> rbackboneghost = r + b1 * POS_BACK - a1 * POS_BACK;
		number rbackghostmod = rbackboneghost.module();
		LR_vector<number> rbackboneghostdir = rbackboneghost / rbackghostmod;
		number cosphi3 = rstackdir * (rbackboneghostdir.cross(a1));

		// old code
		/*LR_vector<number> rbackbone = r + q->pos_back - p->pos_back;
		number rbackmod = rbackbone.module();
		LR_vector<number> rbackbonedir = rbackbone / rbackmod;
		number cosphi3 = rstackdir * (rbackbonedir.cross(a1));*/

		// functions called at their relevant arguments
		number f2 = this->_interaction.f2(rstackmod, CXST_F2);
		/*
		 number f4t1 = this->_interaction.f4(t1, CXST_F4_THETA1) + this->_interaction.f4(2 * PI - t1, CXST_F4_THETA1);
		 number f4t4 = this->_interaction.f4(t4, CXST_F4_THETA4);
		 number f4t5 = this->_interaction.f4(t5, CXST_F4_THETA5) + this->_interaction.f4(PI - t5, CXST_F4_THETA5);
		 number f4t6 = this->_interaction.f4(t6, CXST_F4_THETA6) + this->_interaction.f4(PI - t6, CXST_F4_THETA6);
		 */
		number f4t1 = this->_interaction.query_mesh(cost1,
				this->_interaction.mesh_f4[CXST_F4_THETA1]);
		number f4t4 = this->_interaction.query_mesh(cost4,
				this->_interaction.mesh_f4[CXST_F4_THETA4]);
		number f4t5 = this->_interaction.query_mesh(cost5,
				this->_interaction.mesh_f4[CXST_F4_THETA5])
				+ this->_interaction.query_mesh(-cost5,
						this->_interaction.mesh_f4[CXST_F4_THETA5]);
		number f4t6 = this->_interaction.query_mesh(cost6,
				this->_interaction.mesh_f4[CXST_F4_THETA6])
				+ this->_interaction.query_mesh(-cost6,
						this->_interaction.mesh_f4[CXST_F4_THETA6]);
		number f5cosphi3 = this->_interaction.f5(cosphi3, CXST_F5_PHI3);

		number cxst_energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);
		energy += cxst_energy;

		// again, makes sense (see same line for HB)
		if (cxst_energy != 0.) {
			// derivatives called at the relevant arguments
			number f2D = this->_interaction.f2D(rstackmod, CXST_F2);
			/*
			 number f4t1Dsin = -this->_interaction.f4Dsin(t1, CXST_F4_THETA1) + this->_interaction.f4Dsin(2 * PI - t1, CXST_F4_THETA1);
			 number f4t4Dsin =  this->_interaction.f4Dsin(t4, CXST_F4_THETA4);
			 number f4t5Dsin =  this->_interaction.f4Dsin(t5, CXST_F4_THETA5) - this->_interaction.f4Dsin(PI - t5, CXST_F4_THETA5);
			 number f4t6Dsin = -this->_interaction.f4Dsin(t6, CXST_F4_THETA6) + this->_interaction.f4Dsin(PI - t6, CXST_F4_THETA6);*/
			number f4t1Dsin = this->_interaction.query_meshD(cost1,
					this->_interaction.mesh_f4[CXST_F4_THETA1]);
			number f4t4Dsin = -this->_interaction.query_meshD(cost4,
					this->_interaction.mesh_f4[CXST_F4_THETA4]);
			number f4t5Dsin = -this->_interaction.query_meshD(cost5,
					this->_interaction.mesh_f4[CXST_F4_THETA5])
					+ this->_interaction.query_meshD(-cost5,
							this->_interaction.mesh_f4[CXST_F4_THETA5]);
			number f4t6Dsin = this->_interaction.query_meshD(cost6,
					this->_interaction.mesh_f4[CXST_F4_THETA6])
					- this->_interaction.query_meshD(-cost6,
							this->_interaction.mesh_f4[CXST_F4_THETA6]);

			number f5Dcosphi3 = this->_interaction.f5D(cosphi3, CXST_F5_PHI3);

			//f4t1 = 1.;
			//f4t4 = 1.;
			//f4t5 = 1.;
			//f4t6 = 1.;
			//f5cosphi3 = 1.;

			//f4t1Dsin = 0.;
			//f4t4Dsin = 0.;
			//f4t5Dsin = 0.;
			//f4t6Dsin = 0.;
			//f5Dcosphi3 = 0.;

			/*FILE * yo = fopen("yo.dat", "w");
			 for (int kk = 0; kk<10000; kk++)
			 {
			 fprintf(yo,"%lf %lf %lf\n", 0.0002 * (kk - 5000) * 2 * PI, this->_interaction.f5(0.0002 * (kk-5000) * 2 * PI, CXST_F5_PHI3), this->_interaction.f5D(0.0002 * (kk - 5000) * 2 * PI, CXST_F5_PHI3));
			 }
			 fclose(yo);
			 exit(-1);*/

			// RADIAL PART
			force = -rstackdir * f2D * f4t1 * f4t4 * f4t5
					* f4t6 * SQR(f5cosphi3);

			// THETA1; t1 = LRACOS (-a1 * b1);
			LR_vector<number> dir = a1.cross(b1);
			number torquemod = -f2 * f4t1Dsin * f4t4 * f4t5
					* f4t6 * SQR(f5cosphi3);

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// TETA4; t4 = LRACOS (a3 * b3);
			dir = a3.cross(b3);
			torquemod = -f2 * f4t1 * f4t4Dsin * f4t5 * f4t6 * SQR(f5cosphi3);

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// THETA5; t5 = LRACOS ( a3 * rstackdir);
			number fact = f2 * f4t1 * f4t4 * f4t5Dsin * f4t6 * SQR(f5cosphi3);
			//force += fact * (a3 - rstackdir * cos(t5)) / rstackmod;
			force += fact * (a3 - rstackdir * cost5) / rstackmod;
			dir = rstackdir.cross(a3);
			torquep -= dir * fact;

			// THETA6; t6 = LRACOS (-b3 * rstackdir);
			fact = f2 * f4t1 * f4t4 * f4t5 * f4t6Dsin * SQR(f5cosphi3);
			//force += (b3 + rstackdir * cos(t6)) / rstackmod * fact;
			force += (b3 + rstackdir * cost6) / rstackmod * fact;
			dir = rstackdir.cross(b3);

			torqueq += -dir * fact;

			// Cosphi3 (qui son dolori...) (meno male che cosphi4 = cosphi3)
			// Definition used:
			// cosphi3 = gamma * (ra3 * a2b1 - ra2 * a3b1)/ rbackmod
			number gamma = POS_STACK - POS_BACK;
			number gammacub = gamma * gamma * gamma;
			number rbackghostmodcub = rbackghostmod * rbackghostmod * rbackghostmod;
			//number a1b1 = a1 * b1;
			number a2b1 = a2 * b1;
			number a3b1 = a3 * b1;
			number ra1 = rstackdir * a1;
			number ra2 = rstackdir * a2;
			number ra3 = rstackdir * a3;
			number rb1 = rstackdir * b1;

			number parentesi = (ra3 * a2b1 - ra2 * a3b1);

			number dcdr = -gamma * parentesi * (gamma * (ra1 - rb1) + rstackmod)
					/ rbackghostmodcub;
			number dcda1b1 = gammacub * parentesi / rbackghostmodcub;
			number dcda2b1 = gamma * ra3 / rbackghostmod;
			number dcda3b1 = -gamma * ra2 / rbackghostmod;
			number dcdra1 = -SQR(gamma) * parentesi * rstackmod / rbackghostmodcub;
			number dcdra2 = -gamma * a3b1 / rbackghostmod;
			number dcdra3 = gamma * a2b1 / rbackghostmod;
			number dcdrb1 = SQR(gamma) * parentesi * rstackmod / rbackghostmodcub;

			number force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * 2 * f5cosphi3
					* f5Dcosphi3;

			force += -force_c
					* (rstackdir * dcdr
							+ ((a1 - rstackdir * ra1) * dcdra1
									+ (a2 - rstackdir * ra2) * dcdra2
									+ (a3 - rstackdir * ra3) * dcdra3
									+ (b1 - rstackdir * rb1) * dcdrb1)
									/ rstackmod);

			torquep += force_c
					* (rstackdir.cross(a1) * dcdra1
							+ rstackdir.cross(a2) * dcdra2
							+ rstackdir.cross(a3) * dcdra3);
			torqueq += force_c * (rstackdir.cross(b1) * dcdrb1);

			LR_vector<number> puretorque = force_c
					* (a1.cross(b1) * dcda1b1 + a2.cross(b1) * dcda2b1
							+ a3.cross(b1) * dcda3b1);
			torquep -= puretorque;
			torqueq += puretorque;

			// final update of forces and torques for CXST
			p->force -= force;
			q->force += force;

			torquep -= p->pos_stack.cross(force);
			torqueq += q->pos_stack.cross(force);
		}
	}

	// Final update of the torques
	// total torques
	p->torque += p->orientationT * torquep;
	q->torque += q->orientationT * torqueq;

	return energy;
}

template<typename number>
void FFS_MD_CPUBackend<number>::_ffs_compute_forces(void) {
	int neigh;
	Particle<number> *p;

	this->_U = this->_U_hydr = (number) 0;
	for (int i = 0; i < this->_N; i++) {
		p = &this->_particles[i];
		this->_U += _particle_particle_bonded_interaction(p);

		p->prepare_list();
		neigh = p->next_neighbour();
		while (neigh != P_VIRTUAL) {
			this->_U += _ffs_particle_particle_interaction(p,
					&this->_particles[neigh]);
			neigh = p->next_neighbour();
		}
	}
}

template<typename number>
void FFS_MD_CPUBackend<number>::sim_step(llint curr_step) {
	get_time(&this->_timer, 0);

	this->_op.reset();

	get_time(&this->_timer, 2);
	this->_first_step(curr_step);
	get_time(&this->_timer, 3);

	get_time(&this->_timer, 6);
	if (this->_are_lists_old == true) {
		this->_update_lists();
		this->_N_updates++;
	}
	get_time(&this->_timer, 7);

	get_time(&this->_timer, 8);
	this->_ffs_compute_forces();
	this->_second_step();
	get_time(&this->_timer, 9);

	get_time(&this->_timer, 10);
	if (this->_thermostat != this->THERMOSTAT_NO
			&& (curr_step % this->_newtonian_steps == 0)) {
		if (this->_thermostat == this->THERMOSTAT_JOHN)
			this->_activate_john_thermostat();
		else
			this->_activate_refresh_thermostat();
	}
	get_time(&this->_timer, 11);

	get_time(&this->_timer, 1);

	_op.fill_distance_parameters<number>(this->_particles, this->_box_side);

	//cout << "I just stepped and bond parameter is " << _op.get_hb_parameter(0) << " and distance is " << _op.get_distance_parameter(0) << endl;
	process_times(&this->_timer);

	if (this->check_stop_conditions()) {
		SimManager::stop = true;
		this->_IO->log(this->_IO->LOG_INFO,
				"Reached stop conditions, stopping in step %d", curr_step);
	}

}

template<typename number>
//void FFS_MD_CPUBackend<number>::init(ifstream &conf_input) {
void FFS_MD_CPUBackend<number>::init(char conf_filename[256]) {
	//MDBackend<number>::init(conf_input);
	MDBackend<number>::init(conf_filename);

	_op.init_from_file(_order_parameters_file, this->_particles, this->_N,
			this->_IO);
	init_ffs_from_file(_ffs_file);
}

/*
 File format:
 {
 action = stop_and
 condition1 = params_a > 1
 condition2 = params_b >= 4
 condition3 = params_c < 4
 }
 */
template<typename number>
void FFS_MD_CPUBackend<number>::init_ffs_from_file(const char *fname) {

	FILE * fin;
	fin = fopen(fname, "r");
	if (fin == NULL)
		this->_IO->die("Cannot open %s", fname);

	input_file input;
	char type_str[512];
	int type;
	loadInput(&input, fin);
	getInputString(&input, "action", type_str, 1);
	type = -1;

	if (strcmp(type_str, "stop_and") == 0) {
		this->_IO->log(this->_IO->LOG_INFO,
				"FFS simulation will stop when all conditions fulfilled");
		type = 0;
	}

	if (strcmp(type_str, "stop_or") == 0) {
		this->_IO->log(
				this->_IO->LOG_INFO,
				"FFS simulation will stop when any of the conditions is fulfilled");
		type = 1;
	}
	if (type != 0 && type != 1) {
		this->_IO->die("Unsupported type of action ");
	} else {
		this->_ffs_type = type;
	}

	char pairname[512];
	char strexpr[512];
	int pairid = 1;
	sprintf(pairname, "condition%d", pairid);
	while (getInputString(&input, pairname, strexpr, 0) == KEY_FOUND) {
		parsed_expression newexpression;
		//cout << "Parsing " << strexpr << endl;
		if (!newexpression.parse_expression(strexpr, &_op))
			this->_IO->die("Failed to parse expression %s", pairname);

		this->_conditions.push_back(newexpression);
		pairid++;
		sprintf(pairname, "condition%d", pairid);
	}

	fclose(fin);
}

template<typename number>
char * FFS_MD_CPUBackend<number>::get_op_state_str(void) {
	int * state = _op.get_hb_states();
	char * aux;
	aux = (char *) _state_str;
	for (int i = 0; i < _op.get_hb_parameters_count(); i++) {
		sprintf(aux, "%2d ", state[i]);
		aux = (char *) _state_str + strlen(_state_str);
	}

	double * dstate = _op.get_distance_states();
	for (int i = 0; i < _op.get_distance_parameters_count(); i++) {
		sprintf(aux, "%5.2f", dstate[i]);
		aux = (char *) _state_str + strlen(_state_str);
	}
	return _state_str;
}

template<typename number>
void FFS_MD_CPUBackend<number>::print_energy(llint curr_step) {
	this->_IO->print_energy(*this, curr_step);
}

template<typename number>
bool FFS_MD_CPUBackend<number>::check_stop_conditions(void) {

	if (this->_ffs_type == 0) { //and option
		for (vector<parsed_expression>::iterator i = this->_conditions.begin();
				i != this->_conditions.end(); i++) {
			if (!(*i).eval_expression(&(this->_op)))
				return false;
		}
		return true;
	} else if (this->_ffs_type == 1) { // OR option
		for (vector<parsed_expression>::iterator i = this->_conditions.begin();
				i != this->_conditions.end(); i++) {
			if ((*i).eval_expression(&(this->_op)))
				return true;
		}
		return false;
	} else
		return true;

}

template class FFS_MD_CPUBackend<float> ;
template class FFS_MD_CPUBackend<double> ;
