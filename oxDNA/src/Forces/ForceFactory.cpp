/*
 * ForceFactory.cpp
 *
 *  Created on: 15/mar/2013
 *      Author: lorenzo
 */

#include "ForceFactory.h"

using namespace std;

ForceFactory::ForceFactory() {

}

ForceFactory::~ForceFactory() {

}

template<typename number>
void ForceFactory::add_force(input_file &inp, BaseParticle<number> **particles, int N, bool is_CUDA, number * box_side_ptr) {
	char type_str[512];
	int type;
	getInputString (&inp, "type", type_str, 1);
	type = -1;
	if (strncmp (type_str, "string", 512) == 0) type = 0;
	if (strncmp (type_str, "twist", 512) == 0) type = 1;
	if (strncmp (type_str, "trap", 512) == 0) type = 2;
	if (strncmp (type_str, "repulsion_plane", 512) == 0) type = 3;
	if (strncmp (type_str, "repulsion_plane_moving", 512) == 0) type = 4;
	if (strncmp (type_str, "mutual_trap", 512) == 0) type = 5;
	if (strcmp (type_str, "lowdim_trap") == 0) type = 6;
	if (type != 0 && type != 1 && type != 2 && type !=3 && type != 4 && type != 5 && type != 6) throw oxDNAException("force type %s not implemented. Aborting", type_str);
	int npart, tmpi;
	if(type == 0 || type == 1 || type == 2 || type == 3 || type == 4 || type == 5 || type == 6) {
		getInputInt(&inp, "particle", &npart, 1);
		if(npart >= N || npart < -1) throw oxDNAException("There is a force exerted on a non-existent particle (%d)", npart);
	}

	char group[512] = "default";
	getInputString(&inp, "group_name", group, 0);
	string group_name(group);

	switch (type) {
	case 0: {
		// constant rate force
		double rate, F0, dirx, diry, dirz;
		char strdir[512];
		getInputDouble(&inp, "F0", &F0, 1);
		getInputDouble(&inp, "rate", &rate, 1);
		getInputString(&inp, "dir", strdir, 1);
		tmpi = sscanf(strdir, "%lf,%lf,%lf", &dirx, &diry, &dirz);
		if (tmpi != 3) throw oxDNAException("could not parse direction in external_forces_file. Dying badly");

		if(npart > -1) {
			OX_LOG(Logger::LOG_INFO, "--> adding constant rate force f = (%g + %g * t) [%g,%g,%g] on particle %i", F0, rate, dirx, diry, dirz, npart);
			// TODO: we need to use 'natural' units (like ps or something) instead of steps
			ConstantRateForce<number> *extF = new ConstantRateForce<number>((number) F0, (number) rate, LR_vector<number>((number)dirx, (number)diry, (number)dirz));
			extF->set_group_name(group_name);
			particles[npart]->add_ext_force(extF);
		}
		else { //npart =  -1 means this affects all particles
			for(int i = 0; i < N; i++) {
				OX_LOG(Logger::LOG_INFO, "--> adding constant rate force f = (%g + %g * t) [%g,%g,%g] on part %i", F0, rate, dirx, diry, dirz, i);
				ConstantRateForce<number> *extF = new ConstantRateForce<number>((number) F0, (number) rate, LR_vector<number>((number)dirx, (number)diry, (number)dirz));
				extF->set_group_name(group_name);
				particles[i]->add_ext_force(extF);
			}
		}
		break;
	}
	case 1: {
		// Rotating traps
		double base, rate, cx, cy, cz, px, py, pz, ax, ay, az, stiff, mx, my, mz;
		char strdir[512];
		getInputDouble(&inp, "stiff", &stiff, 1);
		getInputDouble(&inp, "rate", &rate, 1);
		getInputDouble(&inp, "base", &base, 1);
		getInputString (&inp, "axis", strdir, 1);
		tmpi = sscanf(strdir, "%lf,%lf,%lf", &ax, &ay, &az);
		if (tmpi != 3) throw oxDNAException("could not parse pos0 external_forces_file for a twist trap. Dying badly");
		getInputString (&inp, "center", strdir, 1);
		tmpi = sscanf(strdir, "%lf,%lf,%lf", &cx, &cy, &cz);
		if (tmpi != 3) throw oxDNAException("could not parse pos0 external_forces_file for a twist trap. Dying badly");
		getInputString (&inp, "pos0", strdir, 1);
		tmpi = sscanf(strdir, "%lf,%lf,%lf", &px, &py, &pz);
		if (tmpi != 3) throw oxDNAException("could not parse pos0 external_forces_file for a twist trap. Dying badly");
		OX_LOG(Logger::LOG_INFO, "--> adding rotating trap (twist) on particle %i", npart);

		if (getInputString (&inp, "mask", strdir, 0) == KEY_FOUND) {
			tmpi = sscanf(strdir, "%lf,%lf,%lf", &mx, &my, &mz);
			if (tmpi != 3) throw oxDNAException("could not parse pos0 external_forces_file for a twist trap. Dying badly");
		}
		else {
			mx = my = mz = 1.;
		}
		OX_LOG(Logger::LOG_INFO, "--> adding rotating trap (twist) on particle %i", npart);
		// TODO: we need to use 'natural' units
		// (like ps or something) instead of steps
		ConstantRateTorque<number> *extF = new ConstantRateTorque<number> ((number)stiff, (number) base, (number) rate, LR_vector<number>((number)cx, (number)cy, (number)cz), LR_vector<number>((number) ax, (number) ay, (number) az), LR_vector<number> ((number)px, (number)py, (number)pz), LR_vector<number>((number)mx, (number)my, (number) mz));
		extF->set_group_name(group_name);
		particles[npart]->add_ext_force(extF);
		break;
	}
	case 2: {
		double rate, pos0x, pos0y, pos0z, stiff, dirx, diry, dirz;
		char strdir[512], posdir[512];
		getInputDouble(&inp, "stiff", &stiff, 1);
		getInputDouble(&inp, "rate", &rate, 1);
		getInputString(&inp, "dir", strdir, 1);
		tmpi = sscanf(strdir, "%lf,%lf,%lf", &dirx, &diry, &dirz);
		if (tmpi != 3) throw oxDNAException("could not parse dir in external_forces_file. Dying badly");
		getInputString(&inp, "pos0", posdir, 1);
		tmpi = sscanf(posdir, "%lf,%lf,%lf", &pos0x, &pos0y, &pos0z);
		if (tmpi != 3)
			throw oxDNAException("could not parse pos0 in external_forces_file. Dying badly");
		OX_LOG(Logger::LOG_INFO, "--> adding MovingTrap with stiffnes %lf and pos=[%g,%g,%g] + (%g * t) [%g,%g,%g] on part %i", stiff, pos0x, pos0y, pos0z, rate, dirx, diry, dirz, npart);
		MovingTrap<number> *extF = new MovingTrap<number>((number) stiff, LR_vector<number> ((number)pos0x, (number) pos0y, (number) pos0z), (number) rate, LR_vector<number>((number)dirx, (number)diry, (number)dirz));
		extF->set_group_name(group_name);
		particles[npart]->add_ext_force(extF);
		break;
	}
	case 3: {
		double   stiff, dirx, diry, dirz;
		double position;
		char strdir[512];
		getInputDouble(&inp, "stiff", &stiff, 1);
		getInputDouble(&inp, "position", &position, 1);
		getInputString(&inp, "dir", strdir, 1);
		tmpi = sscanf(strdir, "%lf,%lf,%lf", &dirx, &diry, &dirz);
		if (tmpi != 3)
			throw oxDNAException("could not parse dir in external_forces_file. Dying badly");
		if (tmpi != 3)
			throw oxDNAException("could not parse pos0 in external_forces_file. Dying badly");
		if(npart > -1) {
			OX_LOG(Logger::LOG_INFO, "--> adding RepulsionPlane with stiffnes %lf and pos=[%g *x + %g * y +  %g * z + d = 0 ]  on particle %i", stiff,  dirx, diry, dirz, npart);
			RepulsionPlane<number> *extF = new RepulsionPlane<number>((number) stiff, (number) position, LR_vector<number>((number)dirx, (number)diry, (number)dirz));
			extF->set_group_name(group_name);
			particles[npart]->add_ext_force(extF);
		}
		else { //npart =  -1 means this affects all particles
			for(int i = 0; i < N; i++) {
				OX_LOG(Logger::LOG_INFO, "--> adding RepulsionPlane with stiffnes %lf and pos=[%g *x + %g * y +  %g * z + d = 0 ]  on particle %i", stiff,  dirx, diry, dirz, i);
				RepulsionPlane<number> *extF = new RepulsionPlane<number>((number) stiff, (number) position, LR_vector<number>((number)dirx, (number)diry, (number)dirz));
				extF->set_group_name(group_name);
				particles[i]->add_ext_force(extF);
			}
		}
		break;
	}
	case 4: {
		double stiff, dirx, diry, dirz;
		int ref_particle;
		char strdir[512];
		getInputDouble(&inp, "stiff", &stiff, 1);
		getInputInt(&inp, "ref_particle", &ref_particle, 1);
		getInputString(&inp, "dir", strdir, 1);
		tmpi = sscanf(strdir, "%lf,%lf,%lf", &dirx, &diry, &dirz);
		if (tmpi != 3)
			throw oxDNAException("could not parse dir in external_forces_file. Dying badly");
		if(npart > -1) {
			OX_LOG(Logger::LOG_INFO, "--> adding RepulsionPlaneMoving with stiffnes %lf and pos=[%g *x + %g * y +  %g * z + d = 0 ]  on particle %i", stiff,  dirx, diry, dirz, npart);
			RepulsionPlaneMoving<number> *extF = new RepulsionPlaneMoving<number>((number) stiff, particles[ref_particle], LR_vector<number>((number)dirx, (number)diry, (number)dirz), box_side_ptr);
			extF->set_group_name(group_name);
			particles[npart]->add_ext_force(extF);
		}
		else { //npart =  -1 means this affects all particles
			for(int i = 0; i < N; i++) {
				OX_LOG(Logger::LOG_INFO, "--> adding RepulsionPlaneMoving with stiffnes %lf and pos=[%g *x + %g * y +  %g * z + d = 0 ]  on particle %i", stiff,  dirx, diry, dirz, i);
				RepulsionPlaneMoving<number> *extF = new RepulsionPlaneMoving<number>((number) stiff, particles[ref_particle], LR_vector<number>((number)dirx, (number)diry, (number)dirz), box_side_ptr);
				extF->set_group_name(group_name);
				particles[i]->add_ext_force(extF);
			}
		}
		break;
	}
	case 5: {
		// mutual trap between two particles
		double stiff, r0;
		int ref_particle;
		int PBC = 0;
		getInputDouble(&inp, "stiff", &stiff, 1);
		getInputDouble(&inp, "r0", &r0, 1);
		getInputInt(&inp, "ref_particle", &ref_particle, 1);
		getInputBoolAsInt(&inp, "PBC", &PBC, 0);
		OX_LOG(Logger::LOG_INFO, "--> adding MutualTrap to particle %i with stiff %lf and r0 = %lf with respect to particle %i", npart, stiff, r0, ref_particle);
		MutualTrap<number> *extF;
		extF = new MutualTrap<number>((number) stiff, (number) r0, particles[ref_particle], box_side_ptr, PBC);
		extF->set_group_name(group_name);
		particles[npart]->add_ext_force(extF);
		break;
	}
	case 6: {
		double rate, pos0x, pos0y, pos0z, stiff, dirx, diry, dirz;
		int visX,visY,visZ;
		char strdir[512], posdir[512],strvis[512];
		getInputString(&inp, "visibility", strvis, 1);
		tmpi = sscanf(strvis, "%d,%d,%d",  &visX, &visY, &visZ);
		if (tmpi != 3)
			throw oxDNAException("could not parse dir in external_forces_file. Dying badly");
		getInputDouble(&inp, "stiff", &stiff, 1);
		getInputDouble(&inp, "rate", &rate, 1);
		getInputString(&inp, "dir", strdir, 1);
		tmpi = sscanf(strdir, "%lf,%lf,%lf", &dirx, &diry, &dirz);
		if (tmpi != 3)
			throw oxDNAException("could not parse dir in external_forces_file. Dying badly");
		getInputString(&inp, "pos0", posdir, 1);
		tmpi = sscanf(posdir, "%lf,%lf,%lf", &pos0x, &pos0y, &pos0z);
		if (tmpi != 3)
			throw oxDNAException("could not parse dir in external_forces_file. Dying badly");
		OX_LOG(Logger::LOG_INFO, "--> adding LowdimMovingTrap with stiffness %lf and pos=[%g,%g,%g] + (%g * t) [%g,%g,%g] on part %i and visX=%i visY=%i visZ=%i", stiff, pos0x, pos0y, pos0z, rate, dirx, diry, dirz, npart, visX ? 1:0, visY ? 1:0, visZ ? 1:0);

		LowdimMovingTrap<number> *extF = new LowdimMovingTrap<number>((number) stiff, LR_vector<number> ((number)pos0x, (number) pos0y, (number) pos0z), (number) rate, LR_vector<number>((number)dirx, (number)diry, (number)dirz),visX,visY,visZ);
		particles[npart]->add_ext_force(extF);
		break;
	}
	default:
		OX_LOG(Logger::LOG_INFO, "Probably should't reach this point. Hoping for the best");
		break;
	}
}

template void ForceFactory::add_force<float>(input_file &inp, BaseParticle<float> **particles, int N, bool is_CUDA, float * box_side_ptr);
template void ForceFactory::add_force<double>(input_file &inp, BaseParticle<double> **particles, int N, bool is_CUDA, double * box_side_ptr);
