/*
 * ElasticConstantTensor.cpp
 *
 *  Created on: Oct 30, 2013
 *      Author: Lorenzo
 */

#include "MGAnalysis.h"

#include "Utilities/Utils.h"

#define QUICKHULL_IMPLEMENTATION
#include "quickhull.h"

template<typename number>
MGAnalysis<number>::MGAnalysis() :
				_exponent(6),
				_alpha(0.),
				_volume_only(true),
				_rg_only(false) {

}

template<typename number>
MGAnalysis<number>::~MGAnalysis() {

}

template<typename number>
void MGAnalysis<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	string inter;
	getInputString(&sim_inp, "interaction_type", inter, 1);
	if(inter != "MGInteraction") throw oxDNAException("ElasticConstantTensor is not compatible with the interaction '%s'", inter.c_str());

	getInputInt(&sim_inp, "MG_n", &_exponent, 1);

	number rcut;
	getInputNumber(&sim_inp, "MG_rcut", &rcut, 1);
	_sqr_rcut = SQR(rcut);

	getInputNumber(&sim_inp, "MG_alpha", &_alpha, 0);

	number rfene;
	getInputNumber(&sim_inp, "MG_rfene", &rfene, 1);
	_sqr_rfene = SQR(rfene);

	getInputNumber(&sim_inp, "T", &_T, 1);

	getInputBool(&my_inp, "volume_only", &_volume_only, 0);
	getInputBool(&my_inp, "rg_only", &_rg_only, 0);

	if(_volume_only && _rg_only) {
		throw oxDNAException("MGAnalysis: volume_only and rg_only are incompatible");
	}
}

template<typename number>
void MGAnalysis<number>::init(ConfigInfo<number> &config_info) {
	number rep_rcut = pow(2., 1. / _exponent);
	_sqr_rep_rcut = SQR(rep_rcut);

	_gamma = M_PI / (2.25 - pow(2., 1. / 3.));
	_beta = 2 * M_PI - 2.25 * _gamma;
}

template<typename number>
pair<number, number> MGAnalysis<number>::_lame_coefficients() {
	vector<ParticlePair<number> > pairs = this->_config_info.lists->get_potential_interactions();

	number lambda = 0.;
	number mu = 0.;
	typename vector<ParticlePair<number> >::iterator it;
	for(it = pairs.begin(); it != pairs.end(); it++) {
		BaseParticle<number> *p = (*it).first;
		BaseParticle<number> *q = (*it).second;
		LR_vector<number> r = this->_config_info.box->min_image(p->pos, q->pos);
		number r_sqr = r.norm();
		number r_mod = sqrt(r_sqr);

		number first_der = 0.;
		number second_der = 0.;

		if(p->is_bonded(q)) {
			first_der += 30 * _sqr_rfene * r_mod / (_sqr_rfene - r_sqr);
			second_der += 30 * _sqr_rfene / (_sqr_rfene - r_sqr) * (1. + 2 * r_sqr / (_sqr_rfene - r_sqr));
		}

		if(r_sqr < _sqr_rcut) {
			if(r_sqr < _sqr_rep_rcut) {
				first_der += 4. * (6. * pow(r_mod, -7.) - 12. * pow(r_mod, -13.));
				second_der += 4. * (12. * 13. * pow(r_mod, -14.) - 6. * 7. * pow(r_mod, -8.));
			}
			else {
				number sin_part = sin(_gamma * r_sqr + _beta);
				number cos_part = cos(_gamma * r_sqr + _beta);
				first_der += -_alpha * _gamma * r_mod * sin_part;
				second_der += -(_alpha * _gamma * sin_part + 2 * _alpha * SQR(_gamma) * r_sqr * cos_part);
			}
		}

		number factor = second_der - first_der / r_mod;
		for(int d1 = 0; d1 < 3; d1++) {
			for(int d2 = 0; d2 < 3; d2++) {
				if(d1 != d2) {
					number contrib = factor * SQR(r[d1]) * SQR(r[d2]) / r_sqr;
					lambda += contrib;
					mu += contrib;
				}
			}
		}
	}

	lambda /= 6.;
	mu /= 6.;
	mu += 2 * _T * *this->_config_info.N;

	return pair<number, number>(lambda, mu);
}

template<typename number>
number MGAnalysis<number>::_volume() {
	LR_vector<number> com = _com();
	int N = *this->_config_info.N;

	vector<qh_vertex_t> vertices(N);
	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		LR_vector<number> p_pos = this->_config_info.box->get_abs_pos(p) - com;

		vertices[i].x = p_pos.x;
		vertices[i].y = p_pos.y;
		vertices[i].z = p_pos.z;
	}

	qh_mesh_t mesh = qh_quickhull3d(vertices.data(), N);

	number volume = 0.;
	for(int i = 0, j = 0; i < (int)mesh.nindices; i += 3, j++) {
		LR_vector<number> p1(mesh.vertices[mesh.indices[i + 0]].x, mesh.vertices[mesh.indices[i + 0]].y, mesh.vertices[mesh.indices[i + 0]].z);
		LR_vector<number> p2(mesh.vertices[mesh.indices[i + 1]].x, mesh.vertices[mesh.indices[i + 1]].y, mesh.vertices[mesh.indices[i + 1]].z);
		LR_vector<number> p3(mesh.vertices[mesh.indices[i + 2]].x, mesh.vertices[mesh.indices[i + 2]].y, mesh.vertices[mesh.indices[i + 2]].z);

		volume += (p1 * (p2.cross(p3))) / 6.;
	}

	qh_free_mesh(mesh);

	return volume;
}

void eigen_decomposition(double A[3][3], double V[3][3], double d[3]);

template<typename number>
LR_vector<number> MGAnalysis<number>::_com() {
	LR_vector<number> com;
	int N = *this->_config_info.N;
	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		com += this->_config_info.box->get_abs_pos(p);
	}
	return com / N;
}

template<typename number>
vector<number> MGAnalysis<number>::_rg_eigenvalues() {
	vector<number> res;
	LR_vector<number> com = _com();
	int N = *this->_config_info.N;

	double IM[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		LR_vector<number> i_pos = this->_config_info.box->get_abs_pos(p) - com;
		
		IM[0][0] += SQR(i_pos[1]) + SQR(i_pos[2]);
		IM[0][1] += -i_pos[0] * i_pos[1];
		IM[0][2] += -i_pos[0] * i_pos[2];
		
		IM[1][1] += SQR(i_pos[0]) + SQR(i_pos[2]);
		IM[1][2] += -i_pos[1] * i_pos[2];
		
		IM[2][2] += SQR(i_pos[0]) + SQR(i_pos[1]);
	}
	IM[1][0] = IM[0][1];
	IM[2][0] = IM[0][2];
	IM[2][1] = IM[1][2];
	for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++) IM[i][j] /= N;

	double EV[3][3];
	double val[3];
	eigen_decomposition(IM, EV, val);

	res.push_back(sqrt(val[0]));
	res.push_back(sqrt(val[1]));
	res.push_back(sqrt(val[2]));

	LR_vector<number> EVs[3] = {
		LR_vector<number>(EV[0][0], EV[0][1], EV[0][2]),
		LR_vector<number>(EV[1][0], EV[1][1], EV[1][2]),
		LR_vector<number>(EV[2][0], EV[2][1], EV[2][2])
	};
	EVs[0].normalize();
	EVs[1].normalize();
	EVs[2].normalize();

	LR_vector<number> max_along_EVs(-1.e6, -1.e6, -1.e6);
	LR_vector<number> min_along_EVs(1.e6, 1.e6, 1.e6);

	for(int i = 0; i < *this->_config_info.N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		LR_vector<number> p_pos = this->_config_info.box->get_abs_pos(p) - com;

		for(int d = 0; d < 3; d++) {
			number abs = p_pos*EVs[d];
			if(abs > max_along_EVs[d]) max_along_EVs[d] = abs;
			else if(abs < min_along_EVs[d]) min_along_EVs[d] = abs;
		}
	}

	res.push_back(max_along_EVs[0] - min_along_EVs[0]);
	res.push_back(max_along_EVs[1] - min_along_EVs[1]);
	res.push_back(max_along_EVs[2] - min_along_EVs[2]);
	
	return res;
}

template<typename number>
std::string MGAnalysis<number>::get_output_string(llint curr_step) {
	string to_ret;

	if(_volume_only) {
		double volume = _volume();
		to_ret = Utils::sformat("%lf", volume);
	}
	else if(_rg_only) {
		vector<number> eigenvalues = _rg_eigenvalues();
		to_ret = Utils::sformat("%lf %lf %lf %lf %lf %lf", eigenvalues[0], eigenvalues[1], eigenvalues[2], eigenvalues[3], eigenvalues[4], eigenvalues[5]);
	}
	else {
		double volume = _volume();
		pair<number, number> lame = _lame_coefficients();
		number lambda = lame.first / volume;
		number mu = lame.second / volume;
		to_ret = Utils::sformat("%lf %lf", lambda, mu);
	}

	return to_ret;
}

template class MGAnalysis<float> ;
template class MGAnalysis<double> ;

// STUFF TO DIAGONALIZE 3X3 MATRICES
#ifdef MAX
#undef MAX
#endif

#define MAX(a, b) ((a)>(b)?(a):(b))

#define n 3

static double hypot2(double x, double y) {
  return sqrt(x*x+y*y);
}

// Symmetric Householder reduction to tridiagonal form.

static void tred2(double V[n][n], double d[n], double e[n]) {

//  This is derived from the Algol procedures tred2 by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

  for (int j = 0; j < n; j++) {
    d[j] = V[n-1][j];
  }

  // Householder reduction to tridiagonal form.

  for (int i = n-1; i > 0; i--) {

    // Scale to avoid under/overflow.

    double scale = 0.0;
    double h = 0.0;
    for (int k = 0; k < i; k++) {
      scale = scale + fabs(d[k]);
    }
    if (scale == 0.0) {
      e[i] = d[i-1];
      for (int j = 0; j < i; j++) {
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
        V[j][i] = 0.0;
      }
    } else {

      // Generate Householder vector.

      for (int k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      double f = d[i-1];
      double g = sqrt(h);
      if (f > 0) {
        g = -g;
      }
      e[i] = scale * g;
      h = h - f * g;
      d[i-1] = f - g;
      for (int j = 0; j < i; j++) {
        e[j] = 0.0;
      }

      // Apply similarity transformation to remaining columns.

      for (int j = 0; j < i; j++) {
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;
        for (int k = j+1; k <= i-1; k++) {
          g += V[k][j] * d[k];
          e[k] += V[k][j] * f;
        }
        e[j] = g;
      }
      f = 0.0;
      for (int j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      double hh = f / (h + h);
      for (int j = 0; j < i; j++) {
        e[j] -= hh * d[j];
      }
      for (int j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (int k = j; k <= i-1; k++) {
          V[k][j] -= (f * e[k] + g * d[k]);
        }
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
      }
    }
    d[i] = h;
  }

  // Accumulate transformations.

  for (int i = 0; i < n-1; i++) {
    V[n-1][i] = V[i][i];
    V[i][i] = 1.0;
    double h = d[i+1];
    if (h != 0.0) {
      for (int k = 0; k <= i; k++) {
        d[k] = V[k][i+1] / h;
      }
      for (int j = 0; j <= i; j++) {
        double g = 0.0;
        for (int k = 0; k <= i; k++) {
          g += V[k][i+1] * V[k][j];
        }
        for (int k = 0; k <= i; k++) {
          V[k][j] -= g * d[k];
        }
      }
    }
    for (int k = 0; k <= i; k++) {
      V[k][i+1] = 0.0;
    }
  }
  for (int j = 0; j < n; j++) {
    d[j] = V[n-1][j];
    V[n-1][j] = 0.0;
  }
  V[n-1][n-1] = 1.0;
  e[0] = 0.0;
}

// Symmetric tridiagonal QL algorithm.

static void tql2(double V[n][n], double d[n], double e[n]) {

//  This is derived from the Algol procedures tql2, by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

  for (int i = 1; i < n; i++) {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;

  double f = 0.0;
  double tst1 = 0.0;
  double eps = pow(2.0,-52.0);
  for (int l = 0; l < n; l++) {

    // Find small subdiagonal element

    tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
    int m = l;
    while (m < n) {
      if (fabs(e[m]) <= eps*tst1) {
        break;
      }
      m++;
    }

    // If m == l, d[l] is an eigenvalue,
    // otherwise, iterate.

    if (m > l) {
      int iter = 0;
      do {
        iter = iter + 1;  // (Could check iteration count here.)

        // Compute implicit shift

        double g = d[l];
        double p = (d[l+1] - g) / (2.0 * e[l]);
        double r = hypot2(p,1.0);
        if (p < 0) {
          r = -r;
        }
        d[l] = e[l] / (p + r);
        d[l+1] = e[l] * (p + r);
        double dl1 = d[l+1];
        double h = g - d[l];
        for (int i = l+2; i < n; i++) {
          d[i] -= h;
        }
        f = f + h;

        // Implicit QL transformation.

        p = d[m];
        double c = 1.0;
        double c2 = c;
        double c3 = c;
        double el1 = e[l+1];
        double s = 0.0;
        double s2 = 0.0;
        for (int i = m-1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = hypot2(p,e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * (c * g + s * d[i]);

          // Accumulate transformation.

          for (int k = 0; k < n; k++) {
            h = V[k][i+1];
            V[k][i+1] = s * V[k][i] + c * h;
            V[k][i] = c * V[k][i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;

        // Check for convergence.

      } while (fabs(e[l]) > eps*tst1);
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }

  // Sort eigenvalues and corresponding vectors.
  /*
  for (int i = 0; i < n-1; i++) {
    int k = i;
    double p = d[i];
    for (int j = i+1; j < n; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (int j = 0; j < n; j++) {
        p = V[j][i];
        V[j][i] = V[j][k];
        V[j][k] = p;
      }
    }
	}*/
}

void eigen_decomposition(double A[n][n], double V[n][n], double d[n]) {
  double e[n];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      V[i][j] = A[i][j];
    }
  }
  tred2(V, d, e);
  tql2(V, d, e);
}
