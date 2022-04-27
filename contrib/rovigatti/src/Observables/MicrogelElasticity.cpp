/*
 * MicrogelElasticity.cpp
 *
 *  Created on: Mar 13, 2018
 *      Author: Lorenzo
 */

#include "MicrogelElasticity.h"

#include "Utilities/Utils.h"

#include <sstream>

#define QUICKHULL_IMPLEMENTATION
#include "quickhull.h"
#include "diagonalise_3x3.h"

MicrogelElasticity::MicrogelElasticity() {
	_crosslinkers_only = false;
}

MicrogelElasticity::~MicrogelElasticity() {

}

void MicrogelElasticity::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	std::string inter;
	getInputString(&sim_inp, "interaction_type", inter, 1);
	if(inter != "MGInteraction") throw oxDNAException("MicrogelElasticity is not compatible with the interaction '%s'", inter.c_str());

	getInputBool(&my_inp, "crosslinkers_only", &_crosslinkers_only, 0);
	if(_crosslinkers_only) {
		std::string topology;
		getInputString(&sim_inp, "topology", topology, 1);
		std::ifstream topology_file(topology.c_str());
		int N;
		topology_file >> N >> _N_crosslinkers;
		topology_file.close();
	}
}

void MicrogelElasticity::init() {

}

LR_vector MicrogelElasticity::_com() {
	LR_vector com;
	int N = (_crosslinkers_only) ? _N_crosslinkers : _config_info->N();
	for(int i = 0; i < N; i++) {
		BaseParticle *p = _config_info->particles()[i];
		com += _config_info->box->get_abs_pos(p);
	}
	return com / N;
}

std::string MicrogelElasticity::get_output_string(llint curr_step) {
	std::stringstream ss;

	// Convex hull
	LR_vector com = _com();
	int N = (_crosslinkers_only) ? _N_crosslinkers : _config_info->N();

	std::vector<qh_vertex_t> vertices(N);
	for(int i = 0; i < N; i++) {
		BaseParticle *p = _config_info->particles()[i];
		LR_vector p_pos = _config_info->box->get_abs_pos(p) - com;

		vertices[i].x = p_pos.x;
		vertices[i].y = p_pos.y;
		vertices[i].z = p_pos.z;
	}

	qh_mesh_t mesh = qh_quickhull3d(vertices.data(), N);

	number volume = 0.;
	LR_vector ch_com;
	for(int i = 0, j = 0; i < (int) mesh.nindices; i += 3, j++) {
		LR_vector p1(mesh.vertices[mesh.indices[i + 0]].x, mesh.vertices[mesh.indices[i + 0]].y, mesh.vertices[mesh.indices[i + 0]].z);
		LR_vector p2(mesh.vertices[mesh.indices[i + 1]].x, mesh.vertices[mesh.indices[i + 1]].y, mesh.vertices[mesh.indices[i + 1]].z);
		LR_vector p3(mesh.vertices[mesh.indices[i + 2]].x, mesh.vertices[mesh.indices[i + 2]].y, mesh.vertices[mesh.indices[i + 2]].z);

		volume += (p1 * (p2.cross(p3))) / 6.;

		ch_com += (p1 + p2 + p3);
	}
	ss << volume;
	ch_com /= mesh.nindices;

	// Gyration tensor
	double gyration_tensor[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
	for(int i = 0, j = 0; i < (int) mesh.nindices; i += 3, j++) {
		LR_vector p1(mesh.vertices[mesh.indices[i + 0]].x, mesh.vertices[mesh.indices[i + 0]].y, mesh.vertices[mesh.indices[i + 0]].z);
		LR_vector p2(mesh.vertices[mesh.indices[i + 1]].x, mesh.vertices[mesh.indices[i + 1]].y, mesh.vertices[mesh.indices[i + 1]].z);
		LR_vector p3(mesh.vertices[mesh.indices[i + 2]].x, mesh.vertices[mesh.indices[i + 2]].y, mesh.vertices[mesh.indices[i + 2]].z);
		LR_vector triangle_com = (p1 + p2 + p3) / 3. - ch_com;

		gyration_tensor[0][0] += SQR(triangle_com[0]);
		gyration_tensor[0][1] += triangle_com[0] * triangle_com[1];
		gyration_tensor[0][2] += triangle_com[0] * triangle_com[2];

		gyration_tensor[1][1] += SQR(triangle_com[1]);
		gyration_tensor[1][2] += triangle_com[1] * triangle_com[2];

		gyration_tensor[2][2] += SQR(triangle_com[2]);
	}
	gyration_tensor[1][0] = gyration_tensor[0][1];
	gyration_tensor[2][0] = gyration_tensor[0][2];
	gyration_tensor[2][1] = gyration_tensor[1][2];
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			gyration_tensor[i][j] /= mesh.nindices / 3;
		}
	}

	double eigenvectors[3][3];
	double eigenvalues[3];
	eigen_decomposition(gyration_tensor, eigenvectors, eigenvalues);

	volume = 4 * M_PI * sqrt(3) * sqrt(eigenvalues[0]) * sqrt(eigenvalues[1]) * sqrt(eigenvalues[2]);
	ss << " " << volume;
	ss << " " << sqrt(eigenvalues[0]);
	ss << " " << sqrt(eigenvalues[1]);
	ss << " " << sqrt(eigenvalues[2]);

	LR_vector EVs[3] = { LR_vector(eigenvectors[0][0], eigenvectors[0][1], eigenvectors[0][2]), LR_vector(eigenvectors[1][0], eigenvectors[1][1], eigenvectors[1][2]), LR_vector(eigenvectors[2][0], eigenvectors[2][1], eigenvectors[2][2]) };
	EVs[0].normalize();
	EVs[1].normalize();
	EVs[2].normalize();

	qh_free_mesh(mesh);

	return ss.str();
}
