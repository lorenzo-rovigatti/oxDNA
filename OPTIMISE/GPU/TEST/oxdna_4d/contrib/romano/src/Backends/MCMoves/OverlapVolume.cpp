/**
 * @file    OverlapVolume.cpp
 * @date    30/apr/2014
 * @author  flavio
 *
 *
 */
#include "OverlapVolume.h"

/// traslation
template<typename number>
OverlapVolume<number>::OverlapVolume (){
	_pos_old = LR_vector<number> (0., 0., 0.);
	_print_every = (llint) 1e5;
	_rmax = (number) 0.f;

	_mesh_r = -1.f;
	_mesh_w = -1.f;

	_suc = std::vector<std::vector<llint> > ();
	_try = std::vector<std::vector<llint> > ();
}

template<typename number>
OverlapVolume<number>::~OverlapVolume () {

}

template<typename number>
void OverlapVolume<number>::init () {
	BaseMove<number>::init();

	int nbinr = (int) (floor(_rmax / _mesh_r) + 0.001);
	int nbinw = (int) (floor( 2. / _mesh_w) + 0.001);
	_mesh_r = _rmax / nbinr;
	_mesh_w = 2.f / nbinw;
	for (int i = 0; i < nbinr; i ++) {
		std::vector<llint> tmpv1;
		std::vector<llint> tmpv2;
		for (int j = 0; j < nbinw; j ++) {
			tmpv1.push_back(0);
			tmpv2.push_back(0);
		}
		_suc.push_back(tmpv1);
		_try.push_back(tmpv2);
	}

	OX_LOG(Logger::LOG_INFO, "(OverlapVolume.cpp) OverlapVolume Move initiated with rmax = %g (mesh, bins: %g %d - %g %d)", _rmax, _mesh_r, nbinr, _mesh_w, nbinw);

}

template<typename number>
void OverlapVolume<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);

	getInputLLInt (&inp, "print_every", &_print_every, 0);

	getInputNumber(&inp, "rmax", &_rmax, 0);

	getInputNumber(&inp, "mesh_r", &_mesh_r, 1);
	getInputNumber(&inp, "mesh_w", &_mesh_w, 1);


}

template<typename number>
void OverlapVolume<number>::apply (llint curr_step) {

	// we increase the attempted count
	this->_attempted += 1;

	// we select the particle to translate
	int pi = 1;

	BaseParticle<number> * p = this->_Info->particles[pi];

	_pos_old = p->pos;
	_orientation_old = p->orientation;
	_orientationT_old = p->orientationT;

	p->pos.x = 2. * (drand48() - 0.5) * _rmax;
	p->pos.y = 2. * (drand48() - 0.5) * _rmax;
	p->pos.z = 2. * (drand48() - 0.5) * _rmax;

	number t = drand48() * M_PI;
	LR_vector<number> axis = Utils::get_random_vector<number>();

	number sintheta = sin(t);
	number costheta = cos(t);
	number olcos = ((number)1.) - costheta;

	number xyo = axis.x * axis.y * olcos;
	number xzo = axis.x * axis.z * olcos;
	number yzo = axis.y * axis.z * olcos;
	number xsin = axis.x * sintheta;
	number ysin = axis.y * sintheta;
	number zsin = axis.z * sintheta;

	LR_matrix<number> R(axis.x * axis.x * olcos + costheta, xyo - zsin, xzo + ysin,
				xyo + zsin, axis.y * axis.y * olcos + costheta, yzo - xsin,
				xzo - ysin, yzo + xsin, axis.z * axis.z * olcos + costheta);

	p->orientation = p->orientation * R;
	p->orientationT = p->orientation.get_transpose();
	p->set_positions();

	if (p->is_rigid_body()) {
		this->_Info->lists->single_update(p);
		if(!this->_Info->lists->is_updated()) {
			this->_Info->lists->global_update();
		}
	}

	// compute new energy
	this->particle_energy(p);

	BaseParticle<number> *q = this->_Info->particles[0];
	LR_vector<number> dist = this->_Info->box->min_image(q->pos, p->pos);
	number dmod = sqrt(dist.norm());

	bool should_count = true;
	if (dmod >= _rmax) {
		should_count = false;
	}

	if (should_count == true) {
		int rbin = (int) (floor(dmod / _mesh_r) + 0.001);
		int wbin = (int) (floor((1. + p->orientation.v3 * q->orientation.v3) / _mesh_w) + 0.001);

		if (rbin * _mesh_r > _rmax) throw oxDNAException("rbin mesh_r rmax %d %g %g (%g)", rbin, _mesh_r, _rmax,sqrt(this->_Info->box->sqr_min_image_distance(q->pos,p->pos)) );
		if (wbin * _mesh_w > 2.) throw oxDNAException("wbin mesh_w rmax %g %g", wbin, _mesh_w);

		if (rbin < 0) throw oxDNAException ("rbin < 0");
		if (wbin < 0) throw oxDNAException ("rbin < 0");

		_try[rbin][wbin] ++;

		if (this->_Info->interaction->get_is_infinite() == true) {
			_noverlaps ++;
			this->_accepted ++;
			_suc[rbin][wbin] ++;
		}
	}

	p->pos = _pos_old;
	p->orientation = _orientation_old;
	p->orientationT = _orientationT_old;
	p->set_positions();

	if (p->is_rigid_body()) {
		this->_Info->lists->single_update(p);
		if(!this->_Info->lists->is_updated()) {
			this->_Info->lists->global_update();
		}
	}

	if (this->_attempted % _print_every == 0) {
		//number V_out = _rmax * _rmax * _rmax;
		number V_out = (4./3.) * M_PI * _rmax * _rmax * _rmax;
		number W_out = 8 * M_PI * M_PI;
		double f = _noverlaps / (double) this->_attempted;
		OX_LOG(Logger::LOG_INFO, "OverlapVolume: %g %g %g", f, f * V_out * W_out, f * V_out);

		FILE * out = fopen("ovol.dat", "w");
		for (unsigned int i = 0; i < _suc.size(); i ++) {
			for (unsigned int j = 0; j < _suc[i].size(); j ++) {
				//number g = (_suc[i][j] / (double) _try[i][j]) * V_out;
				number g = (_suc[i][j] / (double) _try[i][j]);
				fprintf(out, "%g %g %g %d %d\n", (i + 0.5) * _mesh_r, ((j + 0.5) * _mesh_w) - 1., g, i, j);
			}
			fprintf(out, "\n");
		}
		fclose(out);
	}

	this->_Info->interaction->set_is_infinite(false);

	return;
}

template class OverlapVolume<float>;
template class OverlapVolume<double>;
