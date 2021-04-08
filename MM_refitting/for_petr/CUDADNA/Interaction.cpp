/*
 * Interaction.cpp
 *
 *  Created on: Jun 7, 2011
 *      Author: rovigatti
 */

#include "Interaction.h"
#include <cfloat>
#include "Utils.h"

template <typename number>
Interaction<number>::Interaction() : _average(true) {
	F1_A[0] = HYDR_A;
	F1_A[1] = STCK_A;

	F1_RC[0] = HYDR_RC;
	F1_RC[1] = STCK_RC;

	F1_R0[0] = HYDR_R0;
	F1_R0[1] = STCK_R0;

	F1_BLOW[0] = HYDR_BLOW;
	F1_BLOW[1] = STCK_BLOW;

	F1_BHIGH[0] = HYDR_BHIGH;
	F1_BHIGH[1] = STCK_BHIGH;

	F1_RLOW[0] = HYDR_RLOW;
	F1_RLOW[1] = STCK_RLOW;

	F1_RHIGH[0] = HYDR_RHIGH;
	F1_RHIGH[1] = STCK_RHIGH;

	F1_RCLOW[0] = HYDR_RCLOW;
	F1_RCLOW[1] = STCK_RCLOW;

	F1_RCHIGH[0] = HYDR_RCHIGH;
	F1_RCHIGH[1] = STCK_RCHIGH;

	F2_K[0] = CRST_K;
	F2_K[1] = CXST_K;

	F2_RC[0] = CRST_RC;
	F2_RC[1] = CXST_RC;

	F2_R0[0] = CRST_R0;
	F2_R0[1] = CXST_R0;

	F2_BLOW[0] = CRST_BLOW;
	F2_BLOW[1] = CXST_BLOW;

	F2_BHIGH[0] = CRST_BHIGH;
	F2_BHIGH[1] = CXST_BHIGH;

	F2_RLOW[0] = CRST_RLOW;
	F2_RLOW[1] = CXST_RLOW;

	F2_RHIGH[0] = CRST_RHIGH;
	F2_RHIGH[1] = CXST_RHIGH;

	F2_RCLOW[0] = CRST_RCLOW;
	F2_RCLOW[1] = CXST_RCLOW;

	F2_RCHIGH[0] = CRST_RCHIGH;
	F2_RCHIGH[1] = CXST_RCHIGH;

	F4_THETA_A[0] = STCK_THETA4_A;
	F4_THETA_A[1] = STCK_THETA5_A;

	F4_THETA_A[2] = HYDR_THETA1_A;
	F4_THETA_A[3] = HYDR_THETA2_A;
	F4_THETA_A[4] = HYDR_THETA4_A;
	F4_THETA_A[5] = HYDR_THETA7_A;

	F4_THETA_A[6] = CRST_THETA1_A;
	F4_THETA_A[7]= CRST_THETA2_A;
	F4_THETA_A[8] = CRST_THETA4_A;
	F4_THETA_A[9] = CRST_THETA7_A;

	F4_THETA_A[10] = CXST_THETA1_A;
	F4_THETA_A[11] = CXST_THETA4_A;
	F4_THETA_A[12] = CXST_THETA5_A;

	F4_THETA_B[0] = STCK_THETA4_B;
	F4_THETA_B[1] = STCK_THETA5_B;

	F4_THETA_B[2] = HYDR_THETA1_B;
	F4_THETA_B[3] = HYDR_THETA2_B;
	F4_THETA_B[4] = HYDR_THETA4_B;
	F4_THETA_B[5] = HYDR_THETA7_B;

	F4_THETA_B[6] = CRST_THETA1_B;
	F4_THETA_B[7] = CRST_THETA2_B;
	F4_THETA_B[8] = CRST_THETA4_B;
	F4_THETA_B[9] = CRST_THETA7_B;

	F4_THETA_B[10] = CXST_THETA1_B;
	F4_THETA_B[11] = CXST_THETA4_B;
	F4_THETA_B[12] = CXST_THETA5_B;

	F4_THETA_T0[0] = STCK_THETA4_T0;
	F4_THETA_T0[1] = STCK_THETA5_T0;

	F4_THETA_T0[2] = HYDR_THETA1_T0;
	F4_THETA_T0[3] = HYDR_THETA2_T0;
	F4_THETA_T0[4] = HYDR_THETA4_T0;
	F4_THETA_T0[5] = HYDR_THETA7_T0;

	F4_THETA_T0[6] = CRST_THETA1_T0;
	F4_THETA_T0[7] = CRST_THETA2_T0;
	F4_THETA_T0[8] = CRST_THETA4_T0;
	F4_THETA_T0[9] = CRST_THETA7_T0;

	F4_THETA_T0[10] = CXST_THETA1_T0;
	F4_THETA_T0[11] = CXST_THETA4_T0;
	F4_THETA_T0[12] = CXST_THETA5_T0;

	F4_THETA_TS[0] = STCK_THETA4_TS;
	F4_THETA_TS[1] = STCK_THETA5_TS;

	F4_THETA_TS[2] = HYDR_THETA1_TS;
	F4_THETA_TS[3] = HYDR_THETA2_TS;
	F4_THETA_TS[4] = HYDR_THETA4_TS;
	F4_THETA_TS[5] = HYDR_THETA7_TS;

	F4_THETA_TS[6] = CRST_THETA1_TS;
	F4_THETA_TS[7] = CRST_THETA2_TS;
	F4_THETA_TS[8] = CRST_THETA4_TS;
	F4_THETA_TS[9] = CRST_THETA7_TS;

	F4_THETA_TS[10] = CXST_THETA1_TS;
	F4_THETA_TS[11] = CXST_THETA4_TS;
	F4_THETA_TS[12] = CXST_THETA5_TS;

	F4_THETA_TC[0] = STCK_THETA4_TC;
	F4_THETA_TC[1] = STCK_THETA5_TC;

	F4_THETA_TC[2] = HYDR_THETA1_TC;
	F4_THETA_TC[3] = HYDR_THETA2_TC;
	F4_THETA_TC[4] = HYDR_THETA4_TC;
	F4_THETA_TC[5] = HYDR_THETA7_TC;

	F4_THETA_TC[6] = CRST_THETA1_TC;
	F4_THETA_TC[7] = CRST_THETA2_TC;
	F4_THETA_TC[8] = CRST_THETA4_TC;
	F4_THETA_TC[9] = CRST_THETA7_TC;

	F4_THETA_TC[10] = CXST_THETA1_TC;
	F4_THETA_TC[11] = CXST_THETA4_TC;
	F4_THETA_TC[12] = CXST_THETA5_TC;

	F5_PHI_A[0] = STCK_PHI1_A;
	F5_PHI_A[1] = STCK_PHI2_A;
	F5_PHI_A[2] = CXST_PHI3_A;
	F5_PHI_A[3] = CXST_PHI4_A;

	F5_PHI_B[0] = STCK_PHI1_B;
	F5_PHI_B[1] = STCK_PHI2_B;
	F5_PHI_B[2] = CXST_PHI3_B;
	F5_PHI_B[3] = CXST_PHI3_B;

	F5_PHI_XC[0] = STCK_PHI1_XC;
	F5_PHI_XC[1] = STCK_PHI2_XC;
	F5_PHI_XC[2] = CXST_PHI3_XC;
	F5_PHI_XC[3] = CXST_PHI4_XC;

	F5_PHI_XS[0] = STCK_PHI1_XS;
	F5_PHI_XS[1] = STCK_PHI2_XS;
	F5_PHI_XS[2] = CXST_PHI3_XS;
	F5_PHI_XS[3] = CXST_PHI4_XS;

	MESH_F4_POINTS[HYDR_F4_THETA1] = HYDR_T1_MESH_POINTS;
	MESH_F4_POINTS[HYDR_F4_THETA2] = HYDR_T2_MESH_POINTS;
	MESH_F4_POINTS[HYDR_F4_THETA4] = HYDR_T4_MESH_POINTS;
	MESH_F4_POINTS[HYDR_F4_THETA7] = HYDR_T7_MESH_POINTS;

	MESH_F4_POINTS[STCK_F4_THETA4] = STCK_T4_MESH_POINTS;
	MESH_F4_POINTS[STCK_F4_THETA5] = STCK_T5_MESH_POINTS;

	MESH_F4_POINTS[CRST_F4_THETA1] = CRST_T1_MESH_POINTS;
	MESH_F4_POINTS[CRST_F4_THETA2] = CRST_T2_MESH_POINTS;
	MESH_F4_POINTS[CRST_F4_THETA4] = CRST_T4_MESH_POINTS;
	MESH_F4_POINTS[CRST_F4_THETA7] = CRST_T7_MESH_POINTS;

	MESH_F4_POINTS[CXST_F4_THETA1] = CXST_T1_MESH_POINTS;
	MESH_F4_POINTS[CXST_F4_THETA4] = CXST_T4_MESH_POINTS;
	MESH_F4_POINTS[CXST_F4_THETA5] = CXST_T5_MESH_POINTS;

	for (int i = 0; i < 13; i++ ) {
		// the order of the interpolation interval extremes is reversed, due to
		// the cosine being monotonically decreasing with increasing x
		int points = MESH_F4_POINTS[i];
		number upplimit = cos(fmax( 0, F4_THETA_T0[i] - F4_THETA_TC[i]));
		number lowlimit = cos(fmin(PI, F4_THETA_T0[i] + F4_THETA_TC[i]));

		if (i != CXST_F4_THETA1)
		  build_mesh(&Interaction<number>::fakef4, &Interaction<number>::fakef4D, (void *)(&i), points, lowlimit, upplimit, mesh_f4[i]);
		else
		  {
			build_mesh(&Interaction<number>::fakef4sp, &Interaction<number>::fakef4Dsp, (void *)(&i), points, lowlimit, upplimit, mesh_f4[i]);
		  }
		assert(lowlimit < upplimit);
	  }
}

template<typename number>
void Interaction<number>::get_settings(input_file &inp) {
	int avg_seq;
	if(getInputInt(&inp, "use_average_seq", &avg_seq, 0) == KEY_FOUND) {
		_average = (bool) avg_seq;
		if(!_average) getInputString(&inp, "seq_dep_file", _seq_filename, 1);
	}

}

template<typename number>
inline void Interaction<number>::init(number T) {
	// set the default values
	for(int i = 0; i < 5; i++) {
		for(int j = 0; j < 5; j++) {
			// stacking
			F1_EPS[STCK_F1][i][j] = STCK_BASE_EPS + STCK_FACT_EPS * T;
			F1_SHIFT[STCK_F1][i][j] = F1_EPS[STCK_F1][i][j] * SQR(1 - exp(-(STCK_RC - STCK_R0) * STCK_A));

			// HB
			F1_EPS[HYDR_F1][i][j] = HYDR_EPS;
			F1_SHIFT[HYDR_F1][i][j] = F1_EPS[HYDR_F1][i][j] * SQR(1 - exp(-(HYDR_RC - HYDR_R0) * HYDR_A));
		}
	}

	// if we want to use the sequence dependence then we overwrite all the epsilong regarding
	// interactions between true bases (i.e. the interactions between dummy bases or regular
	// bases and dummy bases keeps the default value)
	if(!_average) {
		char key[256];
		float tmp_value, stck_fact_eps;

		input_file seq_file;
		loadInputFile(&seq_file, _seq_filename);
		if(seq_file.state == ERROR) {
			good = false;
			return;
		}

		// stacking
		getInputFloat(&seq_file, "STCK_FACT_EPS", &stck_fact_eps, 1);
		for(int i = 0; i < 4; i++) {
				for(int j = 0; j < 4; j++) {
					sprintf(key, "STCK_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
					getInputFloat(&seq_file, key, &tmp_value, 1);
					F1_EPS[STCK_F1][i][j] =  tmp_value * ( 1.0 - stck_fact_eps + (T*9.0 * stck_fact_eps ) ); //F1_EPS[STCK_F1][i][j] = tmp_value + stck_fact_eps * T;
					F1_SHIFT[STCK_F1][i][j] = F1_EPS[STCK_F1][i][j] * SQR(1 - exp(-(STCK_RC - STCK_R0) * STCK_A));
				}
		}

		// HB
		// if X and Y are two bases which can interact via HB, then the file must containg either HYDR_X_Y
		// or HYDR_Y_X
		if(getInputFloat(&seq_file, "HYDR_A_T", &tmp_value, 0) == KEY_NOT_FOUND) {
			getInputFloat(&seq_file, "HYDR_T_A", &tmp_value, 1);
		}
		F1_EPS[HYDR_F1][N_A][N_T] = F1_EPS[HYDR_F1][N_T][N_A] = tmp_value;
		F1_SHIFT[HYDR_F1][N_A][N_T] = F1_SHIFT[HYDR_F1][N_T][N_A] = F1_EPS[HYDR_F1][N_A][N_T] * SQR(1 - exp(-(HYDR_RC - HYDR_R0) * HYDR_A));

		if(getInputFloat(&seq_file, "HYDR_G_C", &tmp_value, 0) == KEY_NOT_FOUND) {
			getInputFloat(&seq_file, "HYDR_C_G", &tmp_value, 1);
		}
		F1_EPS[HYDR_F1][N_G][N_C] = F1_EPS[HYDR_F1][N_C][N_G] = tmp_value;
		F1_SHIFT[HYDR_F1][N_G][N_C] = F1_SHIFT[HYDR_F1][N_C][N_G] = F1_EPS[HYDR_F1][N_G][N_C] * SQR(1 - exp(-(HYDR_RC - HYDR_R0) * HYDR_A));

		cleanInputFile(&seq_file);
	}

	good = true;
}

template <typename number>
Interaction<number>::~Interaction() {

}

template<typename number>
inline number Interaction<number>::f1(number r, int type, int n3, int n5) {
	number val = (number) 0;
	if(r < F1_RCHIGH[type]) {
		if(r > F1_RHIGH[type]) {
			val = F1_EPS[type][n3][n5] * F1_BHIGH[type] * SQR(r - F1_RCHIGH[type]);
		}
		else if(r > F1_RLOW[type]) {
			number tmp = 1 - exp(-(r - F1_R0[type]) * F1_A[type]);
			val = F1_EPS[type][n3][n5] * SQR(tmp) - F1_SHIFT[type][n3][n5];
		}
		else if(r > F1_RCLOW[type]) {
			val = F1_EPS[type][n3][n5] * F1_BLOW[type] * SQR(r - F1_RCLOW[type]);
		}
	}

	return val;
}

template<typename number>
inline number Interaction<number>::f1D(number r, int type, int n3, int n5) {
	number val = (number) 0;
	if(r < F1_RCHIGH[type]) {
		if(r > F1_RHIGH[type]) {
			val = 2 * F1_BHIGH[type] * (r - F1_RCHIGH[type]);
		}
		else if(r > F1_RLOW[type]) {
			number tmp = exp(-(r - F1_R0[type]) * F1_A[type]);
			val = 2 * (1 - tmp) * tmp * F1_A[type];
		}
		else if(r > F1_RCLOW[type]) {
			val = 2 * F1_BLOW[type] * (r - F1_RCLOW[type]);
		}
	}

	return F1_EPS[type][n3][n5] * val;
}

template<typename number>
inline number Interaction<number>::f2(number r, int type) {
    number val = (number) 0.;
    if (r < F2_RCHIGH[type]) {
	  if (r > F2_RHIGH[type]) {
		  val = F2_K[type] * F2_BHIGH[type] * SQR(r - F2_RCHIGH[type]);
	  }
	  else if (r > F2_RLOW[type]) {
		  val = (F2_K[type] / 2.) * (SQR(r - F2_R0[type]) - SQR(F2_RC[type] - F2_R0[type]));
	  }
	  else if (r > F2_RCLOW[type]) {
		  val = F2_K[type] * F2_BLOW[type] * SQR(r - F2_RCLOW[type]);
	  }
    }
    return val;
}

template<typename number>
inline number Interaction<number>::f2D(number r, int type) {
    number val = (number) 0.;
    if (r < F2_RCHIGH[type]) {
	if (r > F2_RHIGH[type]) {
	    val = 2. * F2_K[type] * F2_BHIGH[type] * (r - F2_RCHIGH[type]);
	}
	else if (r > F2_RLOW[type]) {
	    val = F2_K[type] * (r - F2_R0[type]);
	}
	else if (r > F2_RCLOW[type]) {
	    val = 2. * F2_K[type] * F2_BLOW[type] * (r - F2_RCLOW[type]);
	}
    }
    return val;
}

template<typename number>
inline number Interaction<number>::fakef4sp(number t, void * par) {
	return f4 (acos(t), *((int*)par)) + f4 (2 * PI - acos(t), *((int*)par));
}

template<typename number>
inline number Interaction<number>::fakef4(number t, void * par) {
	return f4 (acos(t), *((int*)par));
}

template<typename number>
inline number Interaction<number>::fakef4Dsp(number t, void * par) {
	return - f4Dsin (acos(t), *((int*)par)) - f4Dsin (2 * PI - acos(t), *((int*)par));
}

template<typename number>
inline number Interaction<number>::fakef4D(number t, void * par) {
	return - f4Dsin (acos(t), *((int*)par));
}

template<typename number>
inline number Interaction<number>::f4(number t, int type) {
	number val = (number) 0;
	t -= F4_THETA_T0[type];
	if(t < 0) t *= -1;

	if(t < F4_THETA_TC[type]) {
		if(t > F4_THETA_TS[type]) {
			// smoothing
			val = F4_THETA_B[type] * SQR(F4_THETA_TC[type] - t);
		}
		else val = (number) 1.f - F4_THETA_A[type] * SQR(t);
	}

	return val;
}

template<typename number>
inline number Interaction<number>::f4D(number t, int type) {
	number val = (number) 0;
	number m = (number) 1;
	t -= F4_THETA_T0[type];
	// this function is a parabola centered in t0. If t < 0 then the value of the function
	// is the same but the value of its derivative has the opposite sign, so m = -1
	if(t < 0) {
		t *= -1;
		m = (number) -1;
	}

	if(t < F4_THETA_TC[type]) {
		if(t > F4_THETA_TS[type]) {
			// smoothing
			val = m * 2 * F4_THETA_B[type] * (t - F4_THETA_TC[type]);
		}
		else val = -m * 2 * F4_THETA_A[type] * t;
	}

	return val;
}

template<typename number>
inline number Interaction<number>::f4Dsin(number t, int type) {
	number val = (number) 0;
	number m = (number) 1;
	number tt0 = t - F4_THETA_T0[type];
	// this function is a parabola centered in t0. If t < 0 then the value of the function
	// is the same but the value of its derivative has the opposite sign, so m = -1
	if(tt0 < 0) {
		tt0 *= -1;
		m = (number) -1;
	}

	if(tt0 < F4_THETA_TC[type]) {
	    	number sint = sin(t);
		if(tt0 > F4_THETA_TS[type]) {
			// smoothing
			val = m * 2 * F4_THETA_B[type] * (tt0 - F4_THETA_TC[type]) / sint;
		}
		else {
		    if(SQR(sint) > 1e-8) val = -m * 2 * F4_THETA_A[type] * tt0 / sint;
		    else val = -m * 2 * F4_THETA_A[type];
		}
	}

	return val;
}

template<typename number>
inline number Interaction<number>::f5(number f, int type) {
	number val = (number) 0;

	if(f > F5_PHI_XC[type]) {
		if(f < F5_PHI_XS[type]) {
			val = F5_PHI_B[type] * SQR(F5_PHI_XC[type] - f);
		}
		else if(f < 0) {
			val = (number) 1.f - F5_PHI_A[type] * SQR(f);
		}
		else val = (number) 1;
	}

	return val;
}

template<typename number>
inline number Interaction<number>::f5D(number f, int type) {
	number val = (number) 0;

	if(f > F5_PHI_XC[type]) {
		if(f < F5_PHI_XS[type]) {
			val = 2 * F5_PHI_B[type] * (f - F5_PHI_XC[type]);
		}
		else if(f < 0) {
			val = (number) -2.f * F5_PHI_A[type] * f;
		}
		else val = (number) 0;
	}

	return val;
}

template<typename number>
number Interaction<number>::_acosd(number x, void *a) {
	if ( x > 0.998) return  (-acos(0.998))/(1. - 0.998);
	return - 1 / ((number) sqrt(1 - x*x));
}

template<typename number>
number Interaction<number>::_acos(number x, void* a) {
	if ( x > 0.998 )
		return acos(0.998) + (- acos(0.998))/(1 - 0.998)*(x - 0.998) ;
	return (number) acos((double)x);
}

template<typename number>
inline number Interaction<number>::query_meshD(number x, Mesh<number> &m) {
	if (x < m.xlow) return (number) m.B[0];
	if (x >= m.xupp) //return (number) m.B[m.N-1];
		x = m.xupp - FLT_EPSILON;
	int i = (int) ((x - m.xlow)/m.delta);
	number dx = x - m.xlow - m.delta * i;
	return (number) (m.B[i] + (2. * dx * m.C[i] + 3 * dx * dx * m.D[i]) / (m.delta * m.delta));
}

template<typename number>
inline number Interaction<number>::query_mesh(number x, Mesh<number> &m) {
	//if (x < m.xlow) return m.B[0];
	if (x <= m.xlow) return m.A[0];
	if (x >= m.xupp) //return m.A[m.N-1];
		x = m.xupp - FLT_EPSILON;
	int i = (int) ((x - m.xlow)/m.delta);
	number dx = x - m.xlow - m.delta * i;
	return (number) (m.A[i] + dx * (m.B[i] + dx * (m.C[i] + dx * m.D[i]) / (m.delta * m.delta)));
}

template<typename number>
void Interaction<number>::build_mesh(number (Interaction<number>::*f)(number, void*), number (Interaction<number>::*der) (number, void*), void * args, int npoints, number xlow, number xupp, Mesh<number> &m) {
    	assert (xlow < xupp);
	int i;
	number x;

	m.init(npoints);

	number dx = (xupp - xlow) / (number) npoints;
	m.delta = dx;
	m.xlow = xlow;
	m.xupp = xupp;

	number fx0, fx1, derx0, derx1;

	for (i=0; i < npoints + 1; i++) {
		x = xlow + i * dx;

		fx0 = (this->*f)(x, args);
		fx1 = (this->*f)(x + dx, args);
		derx0 = (this->*der)(x, args);
		derx1 = (this->*der)(x + dx, args);

		m.A[i] = fx0;
		m.B[i] = derx0;
		m.D[i] = (2*(fx0 - fx1) + (derx0 + derx1)*dx) / dx;
		m.C[i] = (fx1 - fx0 + (-derx0 - m.D[i]) * dx);
    }
}

template class Interaction<float>;
template class Interaction<double>;
