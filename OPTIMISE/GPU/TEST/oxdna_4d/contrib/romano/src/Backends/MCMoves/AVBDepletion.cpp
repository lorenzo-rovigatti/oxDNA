/**
 * @file    AVBDepletion.cpp
 * @date    30/ott/2018
 * @author  flavio
 *
 */
#include "AVBDepletion.h"

#include <random>  // for poissonian distribution

template<typename number>
AVBDepletion<number>::AVBDepletion (){
	_ntries = -1;
	_sigma_dep = 0.5f;
	_tryvolume = -1.f;
	_mu_dep = (number) 1.f;
	_z = (number) 0.f;
	_generator = std::default_random_engine();
	_rmin = (number) -1.f;
	_rmax = (number) -1.f;
	_delta_rot = (number) -1.f;
	_delta_z = (number) -1.f;
	_Vb = (number) -1.f;
	_rejected_n_too_small[0] = 0;
	_rejected_n_too_small[1] = 0;
	_rejected_overlap[0] = 0;
	_rejected_overlap[1] = 0;
	_rejected_depletion[0] = 0;
	_rejected_depletion[1] = 0;
	_rejected_avb_bias[0] = 0;
	_rejected_avb_bias[1] = 0;
	_split_accepted[0] = 0;
	_split_accepted[1] = 0;
	_skin = (number) 0.f;
	_target_Vb = (number) 0.f;
	_histo_dr = 0.;
	_histo_dt = 0.;
	_nretries = 1;
}

template<typename number>
AVBDepletion<number>::~AVBDepletion () {

}

double clength (double r, double d) {
	return r*r*acos(d/r) - d*sqrt(r*r - d*d);
}

double coverlap (double r1, double r2, double d) {
	if (r2 > r1) {
		double tmp = r1;
		r1 = r2;
		r2 = tmp;
	}
	double d1 = (d*d + r1*r1 - r2*r2) / (2.*d);
	double d2 = (d*d - r1*r1 + r2*r2) / (2.*d);
	return clength(r1,d1) + clength(r2,d2);
}

template<typename number>
void AVBDepletion<number>::init () {
	BaseMove<number>::init();

	_delta_z = _sigma_dep / 2.;

	_tryvolume = (10. + 2. * _sigma_dep) * pow((0.5 + _sigma_dep), 2) * M_PI;

	_z = exp(_mu_dep / this->_T);

	_poisson = std::poisson_distribution<int> (_z*_tryvolume);

	_ntries = _z * _tryvolume; // this will be extracted every time, this is just the average to log

	_skin = 0.;

	_rmax = 1.f + 2. * _sigma_dep - _skin;
	_rmin = 1.f;

	_delta_max = 2. * _sigma_dep;
	_tau_max = (1. - cos(M_PI/4.));
	_zeta_max = 10. / 2.;

	_histo_dr = _delta_max / 100.;
	_histo_dt = _tau_max / 100.;
	_histo_dz = _zeta_max / 100.;

	for (int i = 0; i < 100; i++) {
		for (int j = 0; j < 100; j ++) {
			for (int k = 0; k < 100; k ++) {
				_histo[i][j][k] = 0;
			}
		}
	}
	_gain[0] = _gain[1] = _gain[2] = 0;

	_Vb = 4. * M_PI * M_PI * 2. * (1. - cos(_delta_rot));
	_Vb *= 2. * _delta_z;
	_Vb *= M_PI * (_rmax*_rmax - _rmin * _rmin);

	int N = *this->_Info->N;
	int n = 5;
	number V = this->_Info->box->V() * 8 * M_PI * M_PI;
	number bias_in = ((N - n - 1.) * _Vb / ((n + 1) * (V - _Vb)));
	number bias_out = n * (V - _Vb) /((N - n) * _Vb);

	// let's compute an effective temperature
	double ovolume = coverlap (0.5 + _sigma_dep, 0.5 + _sigma_dep, 1.) - 2. * coverlap(0.5+_sigma_dep, 0.5, 1.);
	ovolume *= (10. + 2.*_sigma_dep);
	double kteq = _z*ovolume;

	//_target_Vb = exp(-1. * kteq) * 3. * V /  (N * 0.5);
	//_target_Vb = exp(-2. * kteq) * V / (N * 0.5);
	_target_Vb = exp(-3. * kteq) * 3.* V / (N * 0.5);
	if (this->_restrict_to_type < 0) throw oxDNAException("Depletion move MUST be restricted to a type");

	OX_LOG(Logger::LOG_INFO, "(AVBDepletion.cpp) AVBDepletion move initiated with rmax= %g, rmin=%g, delta_z=%g, delta_rot=%g", _rmax, _rmin, _delta_z, _delta_rot);
	OX_LOG(Logger::LOG_INFO, "(AVBDepletion.cpp)                               Vb=%g, bias_in=%g, bias_out=%g", _Vb, bias_in, bias_out);
	OX_LOG(Logger::LOG_INFO, "(AVBDepletion.cpp)                               <ntries>=%d, sigma_dep=%g, mu_gas=%g, tryvolume=%g", _ntries, _sigma_dep, _mu_dep, _tryvolume);
	OX_LOG(Logger::LOG_INFO, "(AVBDepletion.cpp)                               restrict_to_type=%d, and probability %g", this->_restrict_to_type, this->prob);
	OX_LOG(Logger::LOG_INFO, "(AVBDepletion.cpp)                               max overlap volume=%g, effective T=%g, target_Vb=%g", ovolume, 1./kteq, _target_Vb);
}

template<typename number>
void AVBDepletion<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);

	getInputNumber (&inp, "delta_rot", &_delta_rot, 1);
	//getInputNumber (&inp, "delta_z", &_delta_z, 1);
	getInputNumber (&inp, "sigma_dep", &_sigma_dep, 1);
	getInputNumber (&inp, "mu_dep", &_mu_dep, 1);
	getInputNumber (&inp, "skin", &_skin, 0);
}

int part_log(int n, int k) {
	if (n < k) throw oxDNAException ("never call me like that, you useless programmer");
	if (n == k) return 1;
	if (n == k + 1) return n;
	else return n * part_log(n - 1, k);
}

template<typename number>
number AVBDepletion<number>::fake_vb (number rmax, number thetamax) {
	return 4. * M_PI * M_PI * 2. * (1. - cos(thetamax)) * M_PI * (rmax * rmax - _rmin * _rmin) * 2. * _delta_z;
}

template<typename number>
bool AVBDepletion<number>::are_close(BaseParticle<number> * pt, BaseParticle<number> * p, bool log) {
	LR_vector<number> dist = this->_Info->box->min_image(pt, p);

	number dx = dist * pt->orientation.v1;
	number dy = dist * pt->orientation.v2;
	number dz = dist * pt->orientation.v3;

	//number cosa = pt->orientation.v3 * p->orientation.v3 / (pt->orientation.v3.module() * p	->orientation.v3.module());
	number cosa = pt->orientation.v3 * p->orientation.v3;

	/*
	number r2 = dx*dx + dy*dy;
	number r = sqrt(r2);
	if (_rmin <= r && r < _rmax && cosa > 1./sqrt(2.) && cosa < 1) {
		int ir = int(floor((sqrt(r2) - _rmin) / _histo_dr) + 0.0001);
		int it = int(floor((1. - fabs(cosa)) / _histo_dt) + 0.0001);
		if (ir >= 100 || ir < 0) throw oxDNAException("disaster ir %d, r=%g, _histo_dr=%g", ir, sqrt(r2), _histo_dr);
		if (it >= 100 || it < 0) throw oxDNAException("disaster it %d", it);
		_histo[ir][it] ++;
	}*/

	number r2 = (dx * dx + dy*dy);
	number r = sqrt(r2);
	number t = LRACOS(fabs(cosa));
	if (r > _rmax && r <= _rmax + _histo_dr) _gain[0] ++;
	if (t > _delta_rot && t <= _delta_rot + _histo_dt) _gain[1] ++;
	if (fabs(dz) > _delta_z && fabs(dz) <= _delta_z + _histo_dz) _gain[2] ++;

	{
		number delta = r - _rmin;
		number tau = 1. - fabs(cosa);
		number zeta = fabs(dz);
		if (delta >= 0 && delta < _delta_max) {
			if (tau >= 0 && tau < _tau_max) {
				if (zeta >= 0 && zeta < _zeta_max) {
					int id = int(floor(delta / _histo_dr) + 0.0001);
					int it = int(floor(tau / _histo_dt) + 0.0001);
					int iz = int(floor(zeta / _histo_dz) + 0.0001);
					_histo[id][it][iz] ++;
				}
			}
		}

	}
	/*
	if (_rmin <= r && r < _rmax && cosa > 1./sqrt(2.) && cosa < 1) {
		int ir = int(floor(r - _rmin) / _histo_dr) + 0.0001);
		int it = int(floor((1. - fabs(cosa)) / _histo_dt) + 0.0001);
		if (ir >= 100 || ir < 0) throw oxDNAException("disaster ir %d, r=%g, _histo_dr=%g", ir, r, _histo_dr);
		if (it >= 100 || it < 0) throw oxDNAException("disaster it %d", it);
		_gain[ir] ++;
	}*/

	//printf ("BB %g %g %g %g (%g)\n", dx, dy, dz, sqrt(r2), acos(cosa));
	if (_rmin * _rmin < r2 && r2 <= _rmax * _rmax ) {
		if (fabs(dz) <= _delta_z) {
			if (fabs(cosa) >= cos(_delta_rot)) {
			//if (fabs(cosa) >= _delta_rot) {
					return true;
			}
		}
	}

	return false;
}

template<typename number>
void AVBDepletion<number>::apply (llint curr_step) {

	number length = 10.;

	// we increase the attempted count
	this->_attempted += 1;
	llint freq = 200000;

	number old_delta_rot = _delta_rot;

	if (this->_attempted % freq == 0 && curr_step >= 10) {
		OX_LOG(Logger::LOG_DEBUG, "rejected because of       <n>: %8.6lf %8.6lf", _rejected_n_too_small[0] / (double) freq, _rejected_n_too_small[1] / (double) freq);
		OX_LOG(Logger::LOG_DEBUG, "rejected because of   overlap: %8.6lf %8.6lf", _rejected_overlap[0] / (double) freq, _rejected_overlap[1] / (double) freq);
		OX_LOG(Logger::LOG_DEBUG, "rejected because of depletion: %8.6lf %8.6lf", _rejected_depletion[0] / (double) freq, _rejected_depletion[1] / (double) freq);
		OX_LOG(Logger::LOG_DEBUG, "rejected because of  avb_bias: %8.6lf %8.6lf", _rejected_avb_bias[0] / (double) freq, _rejected_avb_bias[1] / (double) freq);
		OX_LOG(Logger::LOG_DEBUG, "accepted:                      %8.6lf %8.6lf  (%lld %lld)", _split_accepted[0] /(double) this->_attempted, _split_accepted[1] / (double) this->_attempted, _split_accepted[0], _split_accepted[1]);

		// let's try to fix things
		//double rej_n = _rejected_n_too_small[1] / (double) freq;
		//double rej_o = _rejected_overlap[1] / (double) freq;
		//double rej_d = _rejected_depletion[1] / (double) freq;

		/*
		//if (rej_n < rej_d + rej_o) {
		if (rej_n < rej_d) {
			OX_LOG(Logger::LOG_DEBUG, "   rej_n=%g, rej_o=%g rej_d=%g making delta_rot smaller...", rej_n, rej_o, rej_d);
			//_delta_rot /= this->_rej_fact;
			_delta_rot /= 1.005;
		}
		else {
			OX_LOG(Logger::LOG_DEBUG, "   rej_n=%g, rej_o=%g rej_d=%g making delta_rot bigger....", rej_n, rej_o, rej_d);
			//_delta_rot *= this->_acc_fact;
			_delta_rot *= 1.005;
		}
		if (_delta_rot > 0.20) _delta_rot = 0.2;
		if (_delta_rot < 0.01) _delta_rot = 0.01;
		*/
		// this is what we would do

		/*
		FILE *fp;
		fp = fopen ("histo.txt","w");
		if (fp!=NULL)
		{
			for (int ir = 0; ir < 100; ir ++) {
				for (int it = 0; it < 100; it ++) {
					fprintf(fp,"%g %g %lld\n", _rmin + (ir + 0.5) * _histo_dr, 1. - (it + 0.5) * _histo_dt, _histo[ir][it][dz]);
				}
			}
			fclose (fp);
		}*/

		// comput dhisto / dr, dhisto / dt and dhisto / dz

		/*

		number dhdr = (_gain[0] + 0.0001) / _histo_dr;
		number dhdt = (_gain[1] + 0.0001) / _histo_dt;
		number dhdz = (_gain[2] + 0.0001) / _histo_dz;

		OX_LOG(Logger::LOG_DEBUG, "boundaries: %g %g %g", _rmax, _delta_rot, _delta_z);
		OX_LOG(Logger::LOG_DEBUG, "differentials: %g %g %g", _histo_dr, _histo_dt, _histo_dz);
		OX_LOG(Logger::LOG_DEBUG, "derivatives: %g %g %g (%lld %lld %lld)", dhdr, dhdt, dhdz, _gain[0], _gain[1], _gain[2]);

		// find largest
		int il = 0, is = 2;
		if (dhdr > dhdt && dhdr > dhdz) il = 0;
		if (dhdt > dhdr && dhdt > dhdz) il = 1;
		if (dhdz > dhdr && dhdz > dhdt) il = 2;

		if (dhdr < dhdt && dhdr < dhdz) is = 0;
		if (dhdt < dhdr && dhdt < dhdz) is = 1;
		if (dhdz < dhdr && dhdz < dhdt) is = 2;

		if (is == il) throw oxDNAException("disaster in line %d", __LINE__);
		OX_LOG(Logger::LOG_DEBUG, "largest: %d, smallest: %d", il, is);

		// we increase the largest and decrease the smallest according to _target_Vb;
		if (il == 0) {
			_rmax = _rmax + (_rmax - _rmin)* 1.01;
			if (is == 1) _delta_rot = acos(1. - _target_Vb / (16. * pow(M_PI, 3) * (_rmax*_rmax - _rmin * _rmin) * _delta_z));
			if (is == 2) _delta_z = _target_Vb / (16.*pow(M_PI,3)*(1. - cos(_delta_rot))*(_rmax*_rmax - _rmin * _rmin));
		}
		if (il == 1) {
			_delta_rot = _delta_rot * 1.01;
			if (is == 0) _rmax = sqrt(_target_Vb / (16.*pow(M_PI,3)*(1.-cos(_delta_rot))*_delta_z));
			if (is == 2) {
				//printf("here %g %g %g %g\n", _target_Vb, _delta_rot, _rmax, _rmin);
				_delta_z = _target_Vb / (16.*pow(M_PI,3)*(1. - cos(_delta_rot))*(_rmax*_rmax - _rmin * _rmin));
			}
		}

		if (il == 2) {
			_delta_z = _delta_z * 1.01;
			if (is == 0) _rmax = sqrt(_target_Vb / (16.*pow(M_PI,3)*(1.-cos(_delta_rot))*_delta_z));
			if (is == 1) _delta_rot = acos(1. - _target_Vb / (16. * pow(M_PI, 3) * (_rmax*_rmax - _rmin * _rmin) * _delta_z));
		}

		if (_delta_rot < 1.e-4) _delta_rot = 1.e-4;
		if (_delta_rot > M_PI / 8.) _delta_rot = M_PI / 8.;
		if (_rmax > 1. + 2. * _sigma_dep) _rmax = 1. + 2. * _sigma_dep;
		if (_rmax < 1.001) _rmax = 1.001;
		if (_delta_z < 0.001) _delta_z = 0.001;
		if (_delta_z > length / 2.) _delta_z = length / 2.;

		if (_gain[0] < 500) _histo_dr *= 1.1;
		if (_gain[0] > 1500) _histo_dr *= 0.9;
		if (_gain[1] < 500) _histo_dt *= 1.1;
		if (_gain[1] > 1500) _histo_dt *= 0.9;
		if (_gain[2] < 500) _histo_dr *= 1.1;
		if (_gain[2] > 1500) _histo_dz *= 0.9;

		_gain[0] = _gain[1] = _gain[2] = 0;



		OX_LOG(Logger::LOG_DEBUG, "NEW differentials: %g %g %g", _histo_dr, _histo_dt, _histo_dz);
		OX_LOG(Logger::LOG_DEBUG, "NEW boundaries: %g %g %g", _rmax, _delta_rot, _delta_z);
		*/


		long double total = 0;
		for (int i = 0; i < 100; i ++) {
			for (int j = 0; j < 100; j ++) {
				for (int k = 0; k < 100; k ++) {
					total += _histo[i][j][k];
				}
			}
		}

		/*
		// find best r
		int besti = 0, bestj = 0;
		long double besttotal = -1;
		number bestvb = -1;
		for (int ir = 0; ir < 100; ir ++) {
			int bestji = 0;
			number bestvbi = -1;
			for (int it = 0; it < 100; it ++) {
				number myvb = fake_vb(_rmin + (ir) * _histo_dr, acos(1. - (it) * _histo_dt));
				if (abs(myvb - 10. * _target_Vb) < abs(bestvbi - 10. * _target_Vb)) {
					bestji = it;
					bestvbi = myvb;
				}
			}

			long double mytotal = 0;
			for (int jr = ir; jr < 100; jr ++) {
				for (int jt = bestj; jt < 100; jt ++) {
					mytotal += _histo[jr][jt];
				}
			}
			//printf("\t partial: %d --> %d %g %Lf\n", ir, bestji, bestvbi, mytotal);

			if (mytotal > besttotal && bestvbi > 0.9 * 10. * _target_Vb) {
				besttotal = mytotal;
				bestvb = bestvbi;
				besti = ir;
				bestj = bestji;
			}
		}

		printf ("found best total %Lf (%g %g - %d %d) vb:%g target:%g", besttotal / total, _rmin + besti * _histo_dr, acos(1. - bestj * _histo_dt), besti, bestj, bestvb, _target_Vb);
		_rmax = _rmin + besti * _histo_dr;
		_delta_rot = acos(1. - bestj * _histo_dt);
		*/

		// we go through the histogram and find the best combination
		bool found = false;
		int besti = -1, bestj = -1, bestk = -1;
		long double pmax = -1.;
		for (int i = 0; i < 100; i ++) {
			number delta_max = (i + 0.5) * _histo_dr;
			for (int j = 0; j < 100; j ++) {
				number tau_max = (j + 0.5) * _histo_dt;
				number zeta_max = _target_Vb / (16. * pow(M_PI,3.) * tau_max * (2. * delta_max * _rmin + delta_max * delta_max));
				int k_max = int(floor(zeta_max / _histo_dz) + 0.00001);
				long double myp = 0.;
				for (int ii = 0; ii <= i; ii ++) {
					for (int jj = 0; jj <= j; jj ++) {
						for (int kk = 0; kk <= k_max; kk ++) {
							myp += _histo[ii][jj][kk];
						}
					}
				}
				zeta_max = (k_max + 0.5) * _histo_dz;
				number myvb = 16. * M_PI * M_PI * M_PI * tau_max * (2. * delta_max * _rmin + delta_max * delta_max)	* zeta_max;

				number bias = 3 * this->_Info->box->V() / (*this->_Info->N * myvb);
				number score = (myp / total) * bias;

				if (score < 0.) throw oxDNAException("certo, come no");

				if (score > pmax) {
					pmax = score;
					besti = i;
					bestj = j;
					bestk = k_max;
					found = true;
					//printf ("--> %g\n", myvb);
				}
			}
		}

		OX_LOG(Logger::LOG_DEBUG, "boundaries: %g %g %g, found=%d", _rmax, _delta_rot, _delta_z, found);
		OX_LOG(Logger::LOG_DEBUG, "differentials: %g %g %g", _histo_dr, _histo_dt, _histo_dz);

		if (found) {
			_rmax = _rmin + (besti + 0.5) * _histo_dr;
			_delta_rot = acos(1. - (bestj + 0.5) * _histo_dt);
			_delta_z = (bestk + 0.5) * _histo_dz;
			OX_LOG(Logger::LOG_DEBUG, "NEW boundaries: %g %g %g", _rmax, _delta_rot, _delta_z);

			if (_delta_rot < 1.e-4) _delta_rot = 1.e-4;
			if (_delta_rot > M_PI / 4.) _delta_rot = M_PI / 4.;
			if (_rmax > 1. + 2. * _sigma_dep) _rmax = 1. + 2. * _sigma_dep;
			if (_rmax < 1.001) _rmax = 1.001;
			if (_delta_z < 0.001) _delta_z = 0.001;
			if (_delta_z > length / 2.) _delta_z = length / 2.;
		}
		_Vb = 4. * M_PI * M_PI * 2. * (1. - cos(_delta_rot));
		_Vb *= 2. * _delta_z;
		_Vb *= M_PI * (_rmax * _rmax - _rmin * _rmin);

		/*
		if (_rejected_overlap[1] > _rejected_depletion[1]) {
			_nretries ++;
			OX_LOG(Logger::LOG_INFO, "      increasing _nretries to %d", _nretries);
		}
		//if (_rejected_overlap[1] / (double) freq  < 0.2) {
		else {
			_nretries --;
			OX_LOG(Logger::LOG_INFO, "      decreasing _nretries to %d", _nretries);
		}
		if (_nretries < 1) _nretries = 1;
		if (_nretries > 10) _nretries = 10;
		*/

		int n = 5, N = *this->_Info->N;
		number V = this->_Info->box->V() * 8. * M_PI * M_PI;
		number bias_in = ((N - n - 1.) * _Vb / ((n + 1) * (V - _Vb)));
		number bias_out = n * (V - _Vb) /((N - n) * _Vb);

		_rejected_n_too_small[0] = _rejected_n_too_small[1] = 0;
		_rejected_overlap[0] = _rejected_overlap[1] = 0;
		_rejected_avb_bias[0] = _rejected_avb_bias[1] = 0;
		_rejected_depletion[0] = _rejected_depletion[1] = 0;

		OX_LOG(Logger::LOG_DEBUG, "New vb: %g, bias_in/out=%g/%g; old_delta_rot=%g, new delta_rot=%g", _Vb, bias_in, bias_out, old_delta_rot, _delta_rot);
		//abort();
	}

	// we select the target particle
	int pti = (int) (drand48() * (*this->_Info->N));
	BaseParticle<number> * pt = this->_Info->particles[pti];
	if (this->_restrict_to_type >= 0) {
		while(pt->type != this->_restrict_to_type) {
			pti = (int) (drand48() * (*this->_Info->N));
			pt = this->_Info->particles[pti];
		}
	}
	std::vector<BaseParticle<number> *> target_neighs = this->_Info->lists->get_complete_neigh_list(pt);

	// remove particles that are not in the bonding volume
	typename std::vector<BaseParticle<number> *>::iterator it = target_neighs.begin();
	while (it != target_neighs.end()) {
		if (!are_close(pt, *it, true)) it = target_neighs.erase(it);
		else ++it;
	}

	// we select the particle to move
	int pi = (int) (drand48() * (*this->_Info->N));
	BaseParticle<number> * p = this->_Info->particles[pi];

	bool move_in = drand48() < 0.5;

	if (move_in) {
		/*
		// this turned out to be
		if (target_neighs.size() > 5) {
			_rejected_n_too_small[0] ++;
			return;
		}
		*/
		// choose a NON neighbour of target
		while (std::find(target_neighs.begin(), target_neighs.end(), p) != target_neighs.end() || p == pt) {
			pi = (int) (drand48() * (*this->_Info->N));
			p = this->_Info->particles[pi];
			if (this->_restrict_to_type >= 0) {
				while(p->type != this->_restrict_to_type) {
					pi = (int) (drand48() * (*this->_Info->N));
					p = this->_Info->particles[pi];
				}

			}
			else {
				pi = (int) (drand48() * (*this->_Info->N));
				p = this->_Info->particles[pi];
			}
		}
		if (pi == pti || pt == p) throw oxDNAException("(%d) got same particle %d", __LINE__, p->index);

		//OX_LOG(Logger::LOG_INFO, "moving %d into bonding volume of %d", p->index, pt->index);
	}
	else {
		// choose a neighbour of target
		if (target_neighs.size() == 0) {
			// we may return here, since there is noone to send away
			_rejected_n_too_small[1] ++;
			return;
		}

		if (target_neighs.size() == 1) {
			p = target_neighs[0];
			pi = p->index;
		}
		else {
			pi = (int) (drand48() * (target_neighs.size()));
			p = target_neighs[pi];
			while (p == pt) {
				if (this->_restrict_to_type >= 0) {
					while(p->type != this->_restrict_to_type) {
						pi = (int) (drand48() * (target_neighs.size()));
						p = target_neighs[pi];
					}
				}
				else {
					pi = (int) (drand48() * (target_neighs.size()));
					p = target_neighs[pi];
				}
				//printf ("trying %d...\n", pi);
			}
		}
		if (pt == p) throw oxDNAException("(%d) got same particle pindex %d ptindex %d", __LINE__, p->index, pt->index);
	}

	pi = p->index;
	pti = pt->index;
	if (pi == pti) throw oxDNAException("got same particle...");

	// compute the energy before the move
	number delta_E;
	if (this->_compute_energy_before) delta_E = -this->particle_energy(p);
	else delta_E = (number) 0.f;
	p->set_ext_potential (curr_step, this->_Info->box);
	number delta_E_ext = -p->ext_potential;

	std::vector<BaseParticle<number> *> neighs_old 	= this->_Info->lists->get_complete_neigh_list(p);

	LR_vector<number> pos_old = p->pos;
	LR_matrix<number> orientation_old = p->orientation;
	LR_matrix<number> orientationT_old = p->orientationT;

	number wE = 0.f;

	BaseParticle<number> * psave = new BaseParticle<number> ();

	for (int kt = 0; kt < _nretries; kt ++) {
		this->_Info->interaction->set_is_infinite(false);
		if (move_in) {
			// we move p in the neighbourhood of pt;
			number dpx = -1., dpy = -1., dpz = -1.;
			number cntrl = -1;
			do {
				dpx = (2. * (drand48() - 0.5)) * _rmax;
				dpy = (2. * (drand48() - 0.5)) * _rmax;
				cntrl = dpx * dpx + dpy * dpy;
			} while (!(_rmin * _rmin < cntrl && cntrl < _rmax * _rmax));

			dpz = 2. * (drand48() - 0.5) * _delta_z;

			p->pos = pt->pos + dpx * pt->orientation.v1 + dpy * pt->orientation.v2 + dpz * pt->orientation.v3;

			LR_matrix<number> R;
			number angle;
			if (drand48() < 0.5) angle = acos(cos(_delta_rot) + (1.-cos(_delta_rot))*drand48()); // a [0, delta_rot] --> cos(a) [cos(delta_rot), 1]
			//else angle = acos((-1. + cos(_delta_rot)) * drand48());
			else angle = M_PI - acos(cos(_delta_rot) + (1. - cos(_delta_rot))*drand48());
			//R = Utils::get_random_rotation_matrix_from_angle(angle);

			{
				number phi = drand48() * 2. * M_PI;
				LR_vector<number> axis = cos(phi) * pt->orientation.v1 + sin(phi) * pt->orientation.v2;
				axis.normalize();
				number t = angle;
				number sintheta = sin(t);
				number costheta = cos(t);
				number olcos = 1. - costheta;

				number xyo = axis.x * axis.y * olcos;
				number xzo = axis.x * axis.z * olcos;
				number yzo = axis.y * axis.z * olcos;
				number xsin = axis.x * sintheta;
				number ysin = axis.y * sintheta;
				number zsin = axis.z * sintheta;

				R = LR_matrix<number> (axis.x * axis.x * olcos + costheta, xyo - zsin, xzo + ysin, xyo + zsin, axis.y * axis.y * olcos + costheta, yzo - xsin, xzo - ysin, yzo + xsin, axis.z * axis.z * olcos + costheta);
			}

			p->orientation = pt->orientation * R;
			p->orientationT = p->orientation.get_transpose();
			p->set_positions();

			//printf ("@@ %g %g %g %g (%g - %g)\n", dpx, dpy, dpz, sqrt(cntrl), angle, acos(pt->orientation.v3 * p->orientation.v3 / sqrt(p->orientation.v3.norm()*pt->orientation.v3.norm())));

			//if (are_close (pt, p, false) == false) throw oxDNAException ("investigate...");
			if (are_close(pt, p, false) == false) OX_LOG(Logger::LOG_DEBUG, "they are not close, they should be");
		}
		else {
			// we move p out of the neighbourhood of target
			do {
				LR_vector<number> box_sides = this->_Info->box->box_sides();
				p->pos.x = drand48() * box_sides.x;
				p->pos.y = drand48() * box_sides.y;
				p->pos.z = drand48() * box_sides.z;
				p->orientation = Utils::get_random_rotation_matrix_from_angle<number> (acos(2.*(drand48()-0.5)));
				p->orientationT = p->orientation.get_transpose();
				p->set_positions();
			} while (are_close(pt, p, false) == true);
		}

		this->_Info->lists->single_update(p);
		if(!this->_Info->lists->is_updated()) {
			this->_Info->lists->global_update();
		}

		// chi lo sa
		number newE = this->particle_energy(p);
		if (this->_Info->interaction->get_is_infinite() ==  false) {
			wE = newE;
			psave->pos = p->pos;
			psave->orientation = p->orientation;
			psave->orientationT = p->orientationT;
		}
	}

	p->pos = psave->pos;
	p->orientation = psave->orientation;
	p->orientationT = psave->orientationT;
	delete psave;

	this->_Info->lists->single_update(p);
	if(!this->_Info->lists->is_updated()) {
		this->_Info->lists->global_update();
	}

	// energy after the move
	delta_E += wE;
	//delta_E += this->particle_energy(p);
	p->set_ext_potential(curr_step, this->_Info->box);
	delta_E_ext += p->ext_potential;

	std::vector<BaseParticle<number> *> neighs_new = this->_Info->lists->get_complete_neigh_list(p);

	// if there is no overlap between large particles, we handle the depletant
	//bool depletion_accept = true;
	number depletion_bias = 1.f;
	if (this->_Info->interaction->get_is_infinite() == false) {
		// helper particle pointers for pnew and pold
		BaseParticle<number> * p_new = p;
		BaseParticle<number> * p_old = new BaseParticle<number> ();
		p_old->type = p->type;
		p_old->pos = pos_old;
		p_old->orientation = orientation_old;
		p_old->orientationT = orientationT_old;
		p_old->index = p->index;

		// helper particle pointer for depletant
		BaseParticle<number> * q = new BaseParticle<number> ();
		q->type = this->_restrict_to_type + 1;
		q->index = (*this->_Info->N);

		// for each depletant we test if it overlaps
		//int n_dep = _rho_dep * _tryvolume; // put poissonian here

		_ntries = _poisson(_generator);

		// free volume before the move;
		int ntry = 0, nfv = 0;
		while (ntry < _ntries) {
			number dx = 2.*drand48() - 1.; // between -1 and 1
			number dy = 2.*drand48() - 1.; // between -1 and 1
			number dz = drand48() - 0.5;   // between -0.5 and 0.5;
			while (dx*dx + dy*dy >= 1.) {
				dx = 2.*drand48() - 1.;
				dy = 2.*drand48() - 1.;
			}
			dx = dx * (0.5 + _sigma_dep);
			dy = dy * (0.5 + _sigma_dep);
			dz = dz * (length + 2. * _sigma_dep);

			q->pos = p_new->pos + (p_new->orientation.v1 * dx) +
			                      (p_new->orientation.v2 * dy) +
								  (p_new->orientation.v3 * dz);

			// compare with core
			if (_sigma_dep < 0.5f) {
				if (fabs(dz) > length / 2. - _sigma_dep) {
					if (dx * dx + dy * dy < (0.5 - _sigma_dep) * (0.5 - _sigma_dep)) {
						this->_Info->interaction->set_is_infinite(true); // to skip cycle
					}
				}
			}

			// we check if we overlap with any of the other particles, p_old included
			number junk = this->_Info->interaction->pair_interaction(p_old, q);
			for (it = neighs_new.begin(); it != neighs_new.end() && this->_Info->interaction->get_is_infinite() == false; ++it) {
				junk += this->_Info->interaction->pair_interaction(*it, q);
				if ((*it)->type != this->_restrict_to_type) throw oxDNAException("disaster %d", __LINE__);
			}

			if (this->_Info->interaction->get_is_infinite() == false) nfv ++;
			this->_Info->interaction->set_is_infinite(false);

			ntry ++;
		}

		int n_new = nfv;

		ntry = 0;
		nfv = 0;
		while (ntry < _ntries) {
			number dx = 2.*drand48() - 1.; // between -1 and 1
			number dy = 2.*drand48() - 1.; // between -1 and 1
			number dz = drand48() - 0.5;   // between -0.5 and 0.5;
			while (dx*dx + dy*dy >= 1.) {
				dx = 2.*drand48() - 1.;
				dy = 2.*drand48() - 1.;
			}
			dx = dx * (0.5 + _sigma_dep);
			dy = dy * (0.5 + _sigma_dep);
			dz = dz * (length + 2. * _sigma_dep);

			//q->pos.x = p_old->pos.x + p_old->orientation.v1.x * dx + p_old->orientation.v2.x * dy + p_old->orientation.v3.x * dz;
			//q->pos.y = p_old->pos.y + p_old->orientation.v1.y * dx + p_old->orientation.v2.y * dy + p_old->orientation.v3.y * dz;
			//q->pos.z = p_old->pos.z + p_old->orientation.v1.z * dx + p_old->orientation.v2.z * dy + p_old->orientation.v3.z * dz;

			q->pos = p_old->pos + (p_old->orientation.v1 * dx) +
								  (p_old->orientation.v2 * dy) +
								  (p_old->orientation.v3 * dz);

			// compare with core
  			if (_sigma_dep < 0.5f) {
  				if (fabs(dz) > length / 2. - _sigma_dep) {
  					if (dx * dx + dy * dy < (0.5 - _sigma_dep) * (0.5 - _sigma_dep)) {
  						this->_Info->interaction->set_is_infinite(true); // to skip cycle
  					}
  				}
  			}

			// we check if we overlap with any of the other particles, p_old included
			number junk = this->_Info->interaction->pair_interaction(p_new, q);
			for (it = neighs_old.begin(); it != neighs_old.end() && this->_Info->interaction->get_is_infinite() == false; ++it) {
				junk += this->_Info->interaction->pair_interaction(*it, q);
				if ((*it)->type != this->_restrict_to_type) throw oxDNAException("disaster %d", __LINE__);
			}

			if (this->_Info->interaction->get_is_infinite() == false) nfv ++;
			this->_Info->interaction->set_is_infinite(false);

			ntry ++;
		}

		int n_old = nfv;
		int dn = n_new - n_old;

		number f;
		if (dn >= 0) f = pow(_z * _tryvolume, dn) / part_log(n_old + dn, n_old);
		else f = pow(_z * _tryvolume, dn) * part_log(n_old, n_old + dn);

		depletion_bias = 1./f;

		//printf ("%d %8.6e\n", dn, (double)pow(_z * _tryvolume, dn));
		//if (dn < 0) printf ("HERE %d %g\n\n\n\n", dn, depletion_bias);

		//if (f >= (number) 1.f || f > drand48()) depletion_accept = true;
		//else depletion_accept = false;

		delete q;
		delete p_old;
	}

	int n = target_neighs.size();
	int N = *this->_Info->N;
	number V = this->_Info->box->V() * 8. * M_PI * M_PI;
	number avb_bias;
	if (move_in) {
		avb_bias = ((N - n - 1.) * _Vb)/((n + 1.) * (V - _Vb));
	}
	else {
		avb_bias = (n * (V - _Vb)) / ((N - n) * _Vb);
	}
	// accept or reject?

	//	printf ("%8.6e %8.6e %8.6e\n", delta_E, avb_bias, depletion_bias);
	//if (this->_Info->interaction->get_is_infinite() == false && depletion_bias > 1.) printf ("%8.6e %8.6e %8.6e %d %d\n", delta_E, avb_bias, depletion_bias, move_in, true);

	if (this->_Info->interaction->get_is_infinite() == false &&
		exp(-(delta_E + delta_E_ext) / this->_T) * avb_bias * depletion_bias > drand48() ) {
		// move accepted
        this->_accepted ++;

		if(move_in) _split_accepted[0] ++;
		else _split_accepted[1] ++;

		//printf ("HERE\n");
		//fflush(NULL);

	}
	else {
		p->pos = pos_old;
		p->orientation = orientation_old;
		p->orientationT = orientationT_old;
		p->set_positions();

		if (this->_Info->interaction->get_is_infinite() == true) {
			if (move_in) _rejected_overlap[0] ++;
			else _rejected_overlap[1] ++;
		}
		else {
			if (depletion_bias < avb_bias) {
				if (move_in) _rejected_depletion[0] ++;
				else _rejected_depletion[1] ++;
			}
			else {
				if (move_in) _rejected_avb_bias[0] ++;
				else _rejected_avb_bias[1] ++;
			}
		}

		this->_Info->interaction->set_is_infinite(false);

		this->_Info->lists->single_update(p);
		if(!this->_Info->lists->is_updated()) {
			this->_Info->lists->global_update();
		}
	}

	return;
}

template<typename number>
void AVBDepletion<number>::log_parameters() {
	BaseMove<number>::log_parameters();
	OX_LOG(Logger::LOG_INFO, "AVBDepletion does not change parameters for now...");
}

template class AVBDepletion<float>;
template class AVBDepletion<double>;
