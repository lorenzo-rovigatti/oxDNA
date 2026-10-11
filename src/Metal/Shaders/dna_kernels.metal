/**
 * @file    dna_kernels.metal
 * @brief   Metal kernels for the DNA/DNA2 interaction.
 *
 * This is a faithful port of the CUDA reference implementation in
 *   src/CUDA/Interactions/CUDA_DNA.cuh
 * (functions _excluded_volume, _f1..._f5, _bonded_excluded_volume,
 *  _bonded_part and _particle_particle_DNA_interaction).
 *
 * Forces and torques only: the reported potential energy is recomputed on the
 * CPU by the PotentialEnergy observable, so the kernel does not need to build
 * an energy split. Torque is returned in the particle body frame, exactly like
 * the CUDA kernel (final _vectors_transpose_c_number4_product step), because the
 * Metal velocity-Verlet integrator consumes body-frame angular momenta.
 */

#include "shader_utils.h"

using namespace metal;

// ---------------------------------------------------------------------------
//  Model constants (mirrors src/model.h)
// ---------------------------------------------------------------------------
#define SQR(x)  ((x)*(x))
#define CUB(x)  ((x)*(x)*(x))
#define PI 3.141592653589793f

#define POS_BACK       -0.4f
#define POS_MM_BACK1   -0.3400f
#define POS_MM_BACK2    0.3408f
#define POS_STACK       0.34f
#define POS_BASE        0.4f
#define GAMMA           0.74f

#define FENE_EPS         2.0f
#define FENE_R0_OXDNA    0.7525f
#define FENE_R0_OXDNA2   0.7564f
#define FENE_DELTA       0.25f
#define FENE_DELTA2      0.0625f

#define EXCL_EPS  2.0f
#define EXCL_S1 0.70f
#define EXCL_S2 0.33f
#define EXCL_S3 0.515f
#define EXCL_S4 0.515f
#define EXCL_R1 0.675f
#define EXCL_R2 0.32f
#define EXCL_R3 0.50f
#define EXCL_R4 0.50f
#define EXCL_B1 892.016223343f
#define EXCL_B2 4119.70450017f
#define EXCL_B3 1707.30627298f
#define EXCL_B4 1707.30627298f
#define EXCL_RC1 0.711879214356f
#define EXCL_RC2 0.335388426126f
#define EXCL_RC3 0.52329943261f
#define EXCL_RC4 0.52329943261f

#define HYDR_F1 0
#define STCK_F1 1

#define HYDR_RCLOW  0.276908f
#define HYDR_RCHIGH 0.783775f

#define HYDR_THETA1_A 1.5f
#define HYDR_THETA1_B 4.16038f
#define HYDR_THETA1_T0 0.f
#define HYDR_THETA1_TS 0.7f
#define HYDR_THETA1_TC 0.952381f
#define HYDR_THETA2_A 1.5f
#define HYDR_THETA2_B 4.16038f
#define HYDR_THETA2_T0 0.f
#define HYDR_THETA2_TS 0.7f
#define HYDR_THETA2_TC 0.952381f
#define HYDR_THETA3_A 1.5f
#define HYDR_THETA3_B 4.16038f
#define HYDR_THETA3_T0 0.f
#define HYDR_THETA3_TS 0.7f
#define HYDR_THETA3_TC 0.952381f
#define HYDR_THETA4_A 0.46f
#define HYDR_THETA4_B 0.133855f
#define HYDR_THETA4_T0 PI
#define HYDR_THETA4_TS 0.7f
#define HYDR_THETA4_TC 3.10559f
#define HYDR_THETA7_A 4.f
#define HYDR_THETA7_B 17.0526f
#define HYDR_THETA7_T0 (PI*0.5f)
#define HYDR_THETA7_TS 0.45f
#define HYDR_THETA7_TC 0.555556f
#define HYDR_THETA8_A 4.f
#define HYDR_THETA8_B 17.0526f
#define HYDR_THETA8_T0 (PI*0.5f)
#define HYDR_THETA8_TS 0.45f
#define HYDR_THETA8_TC 0.555556f

#define STCK_THETA4_A 1.3f
#define STCK_THETA4_B 6.4381f
#define STCK_THETA4_T0 0.f
#define STCK_THETA4_TS 0.8f
#define STCK_THETA4_TC 0.961538f
#define STCK_THETA5_A 0.9f
#define STCK_THETA5_B 3.89361f
#define STCK_THETA5_T0 0.f
#define STCK_THETA5_TS 0.95f
#define STCK_THETA5_TC 1.16959f
#define STCK_THETA6_A 0.9f
#define STCK_THETA6_B 3.89361f
#define STCK_THETA6_T0 0.f
#define STCK_THETA6_TS 0.95f
#define STCK_THETA6_TC 1.16959f
#define STCK_F5_PHI1 0
#define STCK_F5_PHI2 1

#define CRST_F2 0
#define CRST_RCLOW  0.45f
#define CRST_RCHIGH 0.7f
#define CRST_THETA1_A 2.25f
#define CRST_THETA1_B 7.00545f
#define CRST_THETA1_T0 (PI - 2.35f)
#define CRST_THETA1_TS 0.58f
#define CRST_THETA1_TC 0.766284f
#define CRST_THETA2_A 1.70f
#define CRST_THETA2_B 6.2469f
#define CRST_THETA2_T0 1.f
#define CRST_THETA2_TS 0.68f
#define CRST_THETA2_TC 0.865052f
#define CRST_THETA3_A 1.70f
#define CRST_THETA3_B 6.2469f
#define CRST_THETA3_T0 1.f
#define CRST_THETA3_TS 0.68f
#define CRST_THETA3_TC 0.865052f
#define CRST_THETA4_A 1.50f
#define CRST_THETA4_B 2.59556f
#define CRST_THETA4_T0 0.f
#define CRST_THETA4_TS 0.65f
#define CRST_THETA4_TC 1.02564f
#define CRST_THETA7_A 1.70f
#define CRST_THETA7_B 6.2469f
#define CRST_THETA7_T0 0.875f
#define CRST_THETA7_TS 0.68f
#define CRST_THETA7_TC 0.865052f
#define CRST_THETA8_A 1.70f
#define CRST_THETA8_B 6.2469f
#define CRST_THETA8_T0 0.875f
#define CRST_THETA8_TS 0.68f
#define CRST_THETA8_TC 0.865052f

#define CXST_F2 1
#define CXST_RCLOW  0.177778f
#define CXST_RCHIGH 0.6222222f
#define CXST_THETA1_A 2.f
#define CXST_THETA1_B 10.9032f
#define CXST_THETA1_T0_OXDNA  (PI - 0.60f)
#define CXST_THETA1_T0_OXDNA2 (PI - 0.25f)
#define CXST_THETA1_TS 0.65f
#define CXST_THETA1_TC 0.769231f
#define CXST_THETA1_SA 20.f
#define CXST_THETA1_SB (PI - 0.1f*(PI - (PI - 0.25f)))
#define CXST_THETA4_A 1.3f
#define CXST_THETA4_B 6.4381f
#define CXST_THETA4_T0 0.f
#define CXST_THETA4_TS 0.8f
#define CXST_THETA4_TC 0.961538f
#define CXST_THETA5_A 0.9f
#define CXST_THETA5_B 3.89361f
#define CXST_THETA5_T0 0.f
#define CXST_THETA5_TS 0.95f
#define CXST_THETA5_TC 1.16959f
#define CXST_THETA6_A 0.9f
#define CXST_THETA6_B 3.89361f
#define CXST_THETA6_T0 0.f
#define CXST_THETA6_TS 0.95f
#define CXST_THETA6_TC 1.16959f
#define CXST_F5_PHI3 2

// ---------------------------------------------------------------------------
//  Runtime parameters (kept in sync with MetalDNAInteraction.mm)
// ---------------------------------------------------------------------------
struct DNAInteractionParams {
    float F1_EPS[50];
    float F1_SHIFT[50];
    float F1_A[2];
    float F1_RC[2];
    float F1_R0[2];
    float F1_BLOW[2];
    float F1_BHIGH[2];
    float F1_RLOW[2];
    float F1_RHIGH[2];
    float F1_RCLOW[2];
    float F1_RCHIGH[2];

    float F2_K[2];
    float F2_RC[2];
    float F2_R0[2];
    float F2_BLOW[2];
    float F2_BHIGH[2];
    float F2_RLOW[2];
    float F2_RHIGH[2];
    float F2_RCLOW[2];
    float F2_RCHIGH[2];

    float F5_PHI_A[4];
    float F5_PHI_B[4];
    float F5_PHI_XC[4];
    float F5_PHI_XS[4];

    float hb_multiplier;
    float T;

    float dh_RC;
    float dh_RHIGH;
    float dh_prefactor;
    float dh_B;
    float dh_minus_kappa;
    int   dh_half_charged_ends;

    int   grooving;
    int   use_oxDNA2_coaxial_stacking;
    int   use_oxDNA2_FENE;
    int   use_debye_huckel;
    int   use_mbf;
    float mbf_xmax;
    float mbf_finf;
    float sqr_rcut;   // COM cutoff^2: beyond this no term contributes
};

struct InitStrandArgs {
    int N;
};

// ---------------------------------------------------------------------------
//  Small helpers
// ---------------------------------------------------------------------------
inline float lracos(float x) {
    return (x >= 1.0f) ? 0.0f : ((x <= -1.0f) ? PI : acos(x));
}

inline float3 stably_normalised(float3 v) {
    float m = fmax(fmax(fabs(v.x), fabs(v.y)), fabs(v.z));
    if(m <= 0.0f) return v;
    float3 res = v / m;
    float res_mod = length(res);
    return (res_mod > 0.0f) ? res / res_mod : res;
}

// quaternion (x,y,z,w) -> body axes, matches CUDA get_vectors_from_quat
inline void get_axes(float4 q, thread float3 &a1, thread float3 &a2, thread float3 &a3) {
    float sqx = q.x * q.x;
    float sqy = q.y * q.y;
    float sqz = q.z * q.z;
    float sqw = q.w * q.w;
    float xy = q.x * q.y;
    float xz = q.x * q.z;
    float xw = q.x * q.w;
    float yz = q.y * q.z;
    float yw = q.y * q.w;
    float zw = q.z * q.w;

    a1 = float3(sqx - sqy - sqz + sqw, 2.0f * (xy + zw),         2.0f * (xz - yw));
    a2 = float3(2.0f * (xy - zw),      -sqx + sqy - sqz + sqw,   2.0f * (yz + xw));
    a3 = float3(2.0f * (xz + yw),      2.0f * (yz - xw),         -sqx - sqy + sqz + sqw);
}

// ---------------------------------------------------------------------------
//  f1..f5 potentials and derivatives (mirrors CUDA_DNA.cuh)
// ---------------------------------------------------------------------------
inline float _f1(float r, int type, int n3, int n5, constant DNAInteractionParams &p) {
    float val = 0.0f;
    if(r < p.F1_RCHIGH[type]) {
        int idx = 25 * type + n3 * 5 + n5;
        if(r > p.F1_RHIGH[type]) {
            val = p.F1_EPS[idx] * p.F1_BHIGH[type] * SQR(r - p.F1_RCHIGH[type]);
        }
        else if(r > p.F1_RLOW[type]) {
            float tmp = 1.0f - exp(-(r - p.F1_R0[type]) * p.F1_A[type]);
            val = p.F1_EPS[idx] * SQR(tmp) - p.F1_SHIFT[idx];
        }
        else if(r > p.F1_RCLOW[type]) {
            val = p.F1_EPS[idx] * p.F1_BLOW[type] * SQR(r - p.F1_RCLOW[type]);
        }
    }
    return val;
}

inline float _f1D(float r, int type, int n3, int n5, constant DNAInteractionParams &p) {
    float val = 0.0f;
    if(r < p.F1_RCHIGH[type]) {
        float eps = p.F1_EPS[25 * type + n3 * 5 + n5];
        if(r > p.F1_RHIGH[type]) {
            val = 2.0f * eps * p.F1_BHIGH[type] * (r - p.F1_RCHIGH[type]);
        }
        else if(r > p.F1_RLOW[type]) {
            float tmp = exp(-(r - p.F1_R0[type]) * p.F1_A[type]);
            val = 2.0f * eps * (1.0f - tmp) * tmp * p.F1_A[type];
        }
        else if(r > p.F1_RCLOW[type]) {
            val = 2.0f * eps * p.F1_BLOW[type] * (r - p.F1_RCLOW[type]);
        }
    }
    return val;
}

inline float _f2(float r, int type, constant DNAInteractionParams &p) {
    float val = 0.0f;
    if(r < p.F2_RCHIGH[type]) {
        if(r > p.F2_RHIGH[type]) {
            val = p.F2_K[type] * p.F2_BHIGH[type] * SQR(r - p.F2_RCHIGH[type]);
        }
        else if(r > p.F2_RLOW[type]) {
            val = (p.F2_K[type] * 0.5f) * (SQR(r - p.F2_R0[type]) - SQR(p.F2_RC[type] - p.F2_R0[type]));
        }
        else if(r > p.F2_RCLOW[type]) {
            val = p.F2_K[type] * p.F2_BLOW[type] * SQR(r - p.F2_RCLOW[type]);
        }
    }
    return val;
}

inline float _f2D(float r, int type, constant DNAInteractionParams &p) {
    float val = 0.0f;
    if(r < p.F2_RCHIGH[type]) {
        if(r > p.F2_RHIGH[type]) {
            val = 2.0f * p.F2_K[type] * p.F2_BHIGH[type] * (r - p.F2_RCHIGH[type]);
        }
        else if(r > p.F2_RLOW[type]) {
            val = p.F2_K[type] * (r - p.F2_R0[type]);
        }
        else if(r > p.F2_RCLOW[type]) {
            val = 2.0f * p.F2_K[type] * p.F2_BLOW[type] * (r - p.F2_RCLOW[type]);
        }
    }
    return val;
}

inline float _f4(float t, float t0, float ts, float tc, float a, float b) {
    float val = 0.0f;
    t = copysign(t - t0, 1.0f);           // |t - t0|
    if(t < tc) {
        val = (t > ts) ? b * SQR(tc - t) : 1.0f - a * SQR(t);
    }
    return val;
}

inline float _f4_pure_harmonic(float t, float a, float b) {
    t -= b;
    return (t < 0.0f) ? 0.0f : a * SQR(t);
}

inline float _f4D(float t, float t0, float ts, float tc, float a, float b) {
    float val = 0.0f;
    t -= t0;
    float m = copysign(1.0f, t);
    t = copysign(t, 1.0f);                 // |t - t0|
    if(t < tc) {
        val = (t > ts) ? 2.0f * m * b * (t - tc) : -2.0f * m * a * t;
    }
    return val;
}

inline float _f4D_pure_harmonic(float t, float a, float b) {
    t -= b;
    return (t < 0.0f) ? 0.0f : 2.0f * a * t;
}

inline float _f5(float f, int type, constant DNAInteractionParams &p) {
    float val = 0.0f;
    if(f > p.F5_PHI_XC[type]) {
        if(f < p.F5_PHI_XS[type]) {
            val = p.F5_PHI_B[type] * SQR(p.F5_PHI_XC[type] - f);
        }
        else if(f < 0.0f) {
            val = 1.0f - p.F5_PHI_A[type] * SQR(f);
        }
        else {
            val = 1.0f;
        }
    }
    return val;
}

inline float _f5D(float f, int type, constant DNAInteractionParams &p) {
    float val = 0.0f;
    if(f > p.F5_PHI_XC[type]) {
        if(f < p.F5_PHI_XS[type]) {
            val = 2.0f * p.F5_PHI_B[type] * (f - p.F5_PHI_XC[type]);
        }
        else if(f < 0.0f) {
            val = -2.0f * p.F5_PHI_A[type] * f;
        }
    }
    return val;
}

// repulsive LJ / smoothed excluded volume, matches CUDA _excluded_volume
inline float3 _excluded_volume(float3 r, float sigma, float rstar, float b, float rc) {
    float rsqr = dot(r, r);
    float3 F = float3(0.0f);
    if(rsqr < SQR(rc)) {
        if(rsqr > SQR(rstar)) {
            float rmod = sqrt(rsqr);
            float rrc = rmod - rc;
            float fmod = 2.0f * EXCL_EPS * b * rrc / rmod;
            F = r * fmod;
        }
        else {
            float lj_part = CUB(SQR(sigma) / rsqr);
            float fmod = 24.0f * EXCL_EPS * (lj_part - 2.0f * SQR(lj_part)) / rsqr;
            F = r * fmod;
        }
    }
    return F;
}

// ---------------------------------------------------------------------------
//  Bonded excluded volume (mirrors CUDA _bonded_excluded_volume<qIsN3>)
// ---------------------------------------------------------------------------
inline void _bonded_excluded_volume(bool qIsN3, float3 r,
                                    float3 n3pos_base, float3 n3pos_back,
                                    float3 n5pos_base, float3 n5pos_back,
                                    thread float3 &F, thread float3 &T) {
    float3 Ftmp;

    Ftmp = _excluded_volume(r + n3pos_base - n5pos_base, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2);
    T += qIsN3 ? cross(n5pos_base, Ftmp) : cross(n3pos_base, Ftmp);
    F += Ftmp;

    Ftmp = _excluded_volume(r + n3pos_back - n5pos_base, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3);
    T += qIsN3 ? cross(n5pos_base, Ftmp) : cross(n3pos_back, Ftmp);
    F += Ftmp;

    Ftmp = _excluded_volume(r + n3pos_base - n5pos_back, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4);
    T += qIsN3 ? cross(n5pos_back, Ftmp) : cross(n3pos_base, Ftmp);
    F += Ftmp;
}

// ---------------------------------------------------------------------------
//  Bonded part: FENE backbone + bonded excluded volume + stacking
//  (mirrors CUDA _bonded_part<qIsN3>)
//    n5* : the 5' partner (parameters n5pos/n5x/n5y/n5z)
//    n3* : the 3' partner (parameters n3pos/n3x/n3y/n3z)
//    r   : n3pos - n5pos (minimum image not needed for bonded neighbours)
// ---------------------------------------------------------------------------
inline void _bonded_part(bool qIsN3, float3 r,
                         int n5type, float3 n5x, float3 n5y, float3 n5z,
                         int n3type, float3 n3x, float3 n3y, float3 n3z,
                         thread float3 &F, thread float3 &T,
                         constant DNAInteractionParams &p) {
    bool grooving = (p.grooving != 0);

    float3 n5pos_back = grooving ? (n5x * POS_MM_BACK1 + n5y * POS_MM_BACK2) : (n5x * POS_BACK);
    float3 n5pos_base = n5x * POS_BASE;
    float3 n5pos_stack = n5x * POS_STACK;

    float3 n3pos_back = grooving ? (n3x * POS_MM_BACK1 + n3y * POS_MM_BACK2) : (n3x * POS_BACK);
    float3 n3pos_base = n3x * POS_BASE;
    float3 n3pos_stack = n3x * POS_STACK;

    float3 rback = r + n3pos_back - n5pos_back;
    float rbackmod = length(rback);
    float rbackr0 = rbackmod - (p.use_oxDNA2_FENE ? FENE_R0_OXDNA2 : FENE_R0_OXDNA);

    float3 Ftmp;
    if(p.use_mbf != 0 && fabs(rbackr0) > p.mbf_xmax) {
        float mbf_fmax = (FENE_EPS * p.mbf_xmax / (FENE_DELTA2 - SQR(p.mbf_xmax)));
        Ftmp = rback * (copysign(1.0f, rbackr0) *
                        ((mbf_fmax - p.mbf_finf) * p.mbf_xmax / fabs(rbackr0) + p.mbf_finf) / rbackmod);
    }
    else {
        Ftmp = rback * ((FENE_EPS * rbackr0 / (FENE_DELTA2 - SQR(rbackr0))) / rbackmod);
    }

    float3 Ttmp = qIsN3 ? cross(n5pos_back, Ftmp) : cross(n3pos_back, Ftmp);

    // EXCLUDED VOLUME (bonded)
    _bonded_excluded_volume(qIsN3, r, n3pos_base, n3pos_back, n5pos_base, n5pos_back, Ftmp, Ttmp);

    if(qIsN3) { F += Ftmp; T += Ttmp; }
    else      { F -= Ftmp; T -= Ttmp; }

    // STACKING
    float3 rstack = r + n3pos_stack - n5pos_stack;
    float rstackmod = length(rstack);
    float3 rstackdir = rstack / rstackmod;

    float3 rbackref = r + n3x * POS_BACK - n5x * POS_BACK;
    float rbackrefmod = length(rbackref);

    float t4 = lracos(dot(n3z, n5z));
    float cost5 = dot(n5z, rstackdir);
    float t5 = lracos(cost5);
    float cost6 = -dot(n3z, rstackdir);
    float t6 = lracos(cost6);
    float cosphi1 = dot(n5y, rbackref) / rbackrefmod;
    float cosphi2 = dot(n3y, rbackref) / rbackrefmod;

    float f1 = _f1(rstackmod, STCK_F1, n3type, n5type, p);
    float f4t4 = _f4(t4, STCK_THETA4_T0, STCK_THETA4_TS, STCK_THETA4_TC, STCK_THETA4_A, STCK_THETA4_B);
    float f4t5 = _f4(PI - t5, STCK_THETA5_T0, STCK_THETA5_TS, STCK_THETA5_TC, STCK_THETA5_A, STCK_THETA5_B);
    float f4t6 = _f4(t6, STCK_THETA6_T0, STCK_THETA6_TS, STCK_THETA6_TC, STCK_THETA6_A, STCK_THETA6_B);
    float f5phi1 = _f5(cosphi1, STCK_F5_PHI1, p);
    float f5phi2 = _f5(cosphi2, STCK_F5_PHI2, p);

    float energy = f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;

    if(energy != 0.0f) {
        float f1D = _f1D(rstackmod, STCK_F1, n3type, n5type, p);
        float f4t4D = _f4D(t4, STCK_THETA4_T0, STCK_THETA4_TS, STCK_THETA4_TC, STCK_THETA4_A, STCK_THETA4_B);
        float f4t5D = _f4D(PI - t5, STCK_THETA5_T0, STCK_THETA5_TS, STCK_THETA5_TC, STCK_THETA5_A, STCK_THETA5_B);
        float f4t6D = _f4D(t6, STCK_THETA6_T0, STCK_THETA6_TS, STCK_THETA6_TC, STCK_THETA6_A, STCK_THETA6_B);
        float f5phi1D = _f5D(cosphi1, STCK_F5_PHI1, p);
        float f5phi2D = _f5D(cosphi2, STCK_F5_PHI2, p);

        // RADIAL
        Ftmp = rstackdir * (energy * f1D / f1);

        // THETA 5
        Ftmp += stably_normalised(n5z - cost5 * rstackdir) * (energy * f4t5D / (f4t5 * rstackmod));

        // THETA 6
        Ftmp += stably_normalised(n3z + cost6 * rstackdir) * (energy * f4t6D / (f4t6 * rstackmod));

        // COS PHI 1
        float ra2 = dot(rstackdir, n5y);
        float ra1 = dot(rstackdir, n5x);
        float rb1 = dot(rstackdir, n3x);
        float a2b1 = dot(n5y, n3x);
        float rbrm3 = SQR(rbackrefmod) * rbackrefmod;

        float dcosphi1dr = (SQR(rstackmod) * ra2 - ra2 * SQR(rbackrefmod)
                            - rstackmod * (a2b1 + ra2 * (-ra1 + rb1)) * GAMMA
                            + a2b1 * (-ra1 + rb1) * SQR(GAMMA)) / rbrm3;
        float dcosphi1dra1 = rstackmod * GAMMA * (rstackmod * ra2 - a2b1 * GAMMA) / rbrm3;
        float dcosphi1dra2 = -rstackmod / rbackrefmod;
        float dcosphi1drb1 = -(rstackmod * GAMMA * (rstackmod * ra2 - a2b1 * GAMMA)) / rbrm3;
        float dcosphi1da1b1 = SQR(GAMMA) * (-rstackmod * ra2 + a2b1 * GAMMA) / rbrm3;
        float dcosphi1da2b1 = GAMMA / rbackrefmod;

        float force_part_phi1 = energy * f5phi1D / f5phi1;

        Ftmp -= (rstackdir * dcosphi1dr
                 + ((n5y - ra2 * rstackdir) * dcosphi1dra2
                    + (n5x - ra1 * rstackdir) * dcosphi1dra1
                    + (n3x - rb1 * rstackdir) * dcosphi1drb1) / rstackmod) * force_part_phi1;

        // COS PHI 2   (p -> b, q -> a)
        ra2 = dot(rstackdir, n3y);
        ra1 = rb1;
        rb1 = dot(rstackdir, n5x);
        a2b1 = dot(n3y, n5x);
        float dcosphi2dr = ((rstackmod * ra2 + a2b1 * GAMMA) * (rstackmod + (rb1 - ra1) * GAMMA)
                            - ra2 * SQR(rbackrefmod)) / rbrm3;
        float dcosphi2dra1 = -rstackmod * GAMMA * (rstackmod * ra2 + a2b1 * GAMMA) / rbrm3;
        float dcosphi2dra2 = -rstackmod / rbackrefmod;
        float dcosphi2drb1 = (rstackmod * GAMMA * (rstackmod * ra2 + a2b1 * GAMMA)) / rbrm3;
        float dcosphi2da1b1 = -SQR(GAMMA) * (rstackmod * ra2 + a2b1 * GAMMA) / rbrm3;
        float dcosphi2da2b1 = -GAMMA / rbackrefmod;

        float force_part_phi2 = energy * f5phi2D / f5phi2;

        Ftmp -= (rstackdir * dcosphi2dr
                 + ((n3y - rstackdir * ra2) * dcosphi2dra2
                    + (n3x - rstackdir * ra1) * dcosphi2dra1
                    + (n5x - rstackdir * rb1) * dcosphi2drb1) / rstackmod) * force_part_phi2;

        Ttmp = qIsN3 ? cross(n5pos_stack, Ftmp) : cross(n3pos_stack, Ftmp);

        // THETA 4
        Ttmp += stably_normalised(cross(n3z, n5z)) * (-energy * f4t4D / f4t4);

        // PHI 1 & PHI 2
        if(qIsN3) {
            Ttmp += (-force_part_phi1 * dcosphi1dra2) * cross(rstackdir, n5y)
                    - cross(rstackdir, n5x) * force_part_phi1 * dcosphi1dra1;
            Ttmp += (-force_part_phi2 * dcosphi2drb1) * cross(rstackdir, n5x);
        }
        else {
            Ttmp += force_part_phi1 * dcosphi1drb1 * cross(rstackdir, n3x);
            Ttmp += force_part_phi2 * dcosphi2dra2 * cross(rstackdir, n3y)
                    + force_part_phi2 * dcosphi2dra1 * cross(rstackdir, n3x);
        }

        Ttmp += force_part_phi1 * dcosphi1da2b1 * cross(n5y, n3x)
                + cross(n5x, n3x) * force_part_phi1 * dcosphi1da1b1;
        Ttmp += force_part_phi2 * dcosphi2da2b1 * cross(n5x, n3y)
                + cross(n5x, n3x) * force_part_phi2 * dcosphi2da1b1;

        if(qIsN3) {
            Ttmp += stably_normalised(cross(rstackdir, n5z)) * (energy * f4t5D / f4t5);
            T += Ttmp;
            F += Ftmp;
        }
        else {
            Ttmp += stably_normalised(cross(rstackdir, n3z)) * (-energy * f4t6D / f4t6);
            T -= Ttmp;
            F -= Ftmp;
        }
    }
}

// ---------------------------------------------------------------------------
//  Non-bonded pair interaction
//  (mirrors CUDA _particle_particle_DNA_interaction)
//   r  : qpos - ppos (minimum image);  a* : p axes;  b* : q axes
//   F,T accumulate onto particle p (IND), returned in the lab frame
// ---------------------------------------------------------------------------
inline void _particle_particle_DNA_interaction(float3 r,
                                               int ptype, float3 a1, float3 a2, float3 a3,
                                               int qtype, float3 b1, float3 b2, float3 b3,
                                               bool p_is_end, bool q_is_end,
                                               thread float3 &F, thread float3 &T,
                                               constant DNAInteractionParams &p) {
    bool grooving = (p.grooving != 0);
    bool use_dh = (p.use_debye_huckel != 0);
    bool use_ox2_cxst = (p.use_oxDNA2_coaxial_stacking != 0);
    int int_type = ptype + qtype;   // WC-complementary pairs sum to 3

    float3 ppos_back = grooving ? (POS_MM_BACK1 * a1 + POS_MM_BACK2 * a2) : (POS_BACK * a1);
    float3 ppos_base = POS_BASE * a1;
    float3 ppos_stack = POS_STACK * a1;

    float3 qpos_back = grooving ? (POS_MM_BACK1 * b1 + POS_MM_BACK2 * b2) : (POS_BACK * b1);
    float3 qpos_base = POS_BASE * b1;
    float3 qpos_stack = POS_STACK * b1;

    // ---- excluded volume (BACK-BACK + the three base/back terms) ----
    float3 Ftmp = float3(0.0f);
    float3 rbackbone = r + qpos_back - ppos_back;
    Ftmp = _excluded_volume(rbackbone, EXCL_S1, EXCL_R1, EXCL_B1, EXCL_RC1);
    float3 Ttmp = cross(ppos_back, Ftmp);
    _bonded_excluded_volume(true, r, qpos_base, qpos_back, ppos_base, ppos_back, Ftmp, Ttmp);
    F += Ftmp;

    // ---- Debye-Huckel ----
    if(use_dh) {
        float rbackmod = length(rbackbone);
        if(rbackmod < p.dh_RC) {
            float3 rbackdir = rbackbone / rbackmod;
            if(rbackmod < p.dh_RHIGH) {
                Ftmp = rbackdir * (-p.dh_prefactor * exp(p.dh_minus_kappa * rbackmod)
                                   * (p.dh_minus_kappa / rbackmod - 1.0f / SQR(rbackmod)));
            }
            else {
                Ftmp = rbackdir * (-2.0f * p.dh_B * (rbackmod - p.dh_RC));
            }
            if(p.dh_half_charged_ends != 0 && p_is_end) Ftmp *= 0.5f;
            if(p.dh_half_charged_ends != 0 && q_is_end) Ftmp *= 0.5f;

            Ttmp -= cross(ppos_back, Ftmp);
            F -= Ftmp;
        }
    }

    // ---- hydrogen bonding ----
    float3 rhydro = r + qpos_base - ppos_base;
    float rhydromodsqr = dot(rhydro, rhydro);
    if(int_type == 3 && SQR(HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(HYDR_RCHIGH)) {
        float hb_multi = 1.0f; // sequence-dependent multiplier handled by F1_EPS
        float rhydromod = sqrt(rhydromodsqr);
        float3 rhydrodir = rhydro / rhydromod;

        float t1 = lracos(-dot(a1, b1));
        float cost2 = -dot(b1, rhydrodir);
        float t2 = lracos(cost2);
        float cost3 = dot(a1, rhydrodir);
        float t3 = lracos(cost3);
        float t4 = lracos(dot(a3, b3));
        float cost7 = -dot(rhydrodir, b3);
        float t7 = lracos(cost7);
        float cost8 = dot(rhydrodir, a3);
        float t8 = lracos(cost8);

        float f1 = hb_multi * _f1(rhydromod, HYDR_F1, ptype, qtype, p);
        float f4t1 = _f4(t1, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC, HYDR_THETA1_A, HYDR_THETA1_B);
        float f4t2 = _f4(t2, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC, HYDR_THETA2_A, HYDR_THETA2_B);
        float f4t3 = _f4(t3, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC, HYDR_THETA3_A, HYDR_THETA3_B);
        float f4t4 = _f4(t4, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC, HYDR_THETA4_A, HYDR_THETA4_B);
        float f4t7 = _f4(t7, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC, HYDR_THETA7_A, HYDR_THETA7_B);
        float f4t8 = _f4(t8, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC, HYDR_THETA8_A, HYDR_THETA8_B);

        float hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

        if(hb_energy < 0.0f) {
            float f1D = hb_multi * _f1D(rhydromod, HYDR_F1, ptype, qtype, p);
            float f4t1D = -_f4D(t1, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC, HYDR_THETA1_A, HYDR_THETA1_B);
            float f4t2D = -_f4D(t2, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC, HYDR_THETA2_A, HYDR_THETA2_B);
            float f4t3D = _f4D(t3, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC, HYDR_THETA3_A, HYDR_THETA3_B);
            float f4t4D = _f4D(t4, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC, HYDR_THETA4_A, HYDR_THETA4_B);
            float f4t7D = -_f4D(t7, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC, HYDR_THETA7_A, HYDR_THETA7_B);
            float f4t8D = _f4D(t8, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC, HYDR_THETA8_A, HYDR_THETA8_B);

            Ftmp = rhydrodir * hb_energy * f1D / f1;

            Ttmp -= stably_normalised(cross(a3, b3)) * (-hb_energy * f4t4D / f4t4);
            Ttmp -= stably_normalised(cross(a1, b1)) * (-hb_energy * f4t1D / f4t1);

            Ftmp -= stably_normalised(b1 + rhydrodir * cost2) * (hb_energy * f4t2D / (f4t2 * rhydromod));

            float part = -hb_energy * f4t3D / f4t3;
            Ftmp -= stably_normalised(a1 - rhydrodir * cost3) * (-part / rhydromod);
            Ttmp += stably_normalised(cross(rhydrodir, a1)) * part;

            Ftmp -= stably_normalised(b3 + rhydrodir * cost7) * (hb_energy * f4t7D / (f4t7 * rhydromod));

            part = -hb_energy * f4t8D / f4t8;
            Ftmp -= stably_normalised(a3 - rhydrodir * cost8) * (-part / rhydromod);
            Ttmp += stably_normalised(cross(rhydrodir, a3)) * part;

            Ttmp += cross(ppos_base, Ftmp);
            F += Ftmp;
        }
    }

    // ---- cross stacking ----
    float3 rcstack = rhydro;
    float rcstackmodsqr = rhydromodsqr;
    if(SQR(CRST_RCLOW) < rcstackmodsqr && rcstackmodsqr < SQR(CRST_RCHIGH)) {
        float rcstackmod = sqrt(rcstackmodsqr);
        float3 rcstackdir = rcstack / rcstackmod;

        float t1 = lracos(-dot(a1, b1));
        float cost2 = -dot(b1, rcstackdir);
        float t2 = lracos(cost2);
        float cost3 = dot(a1, rcstackdir);
        float t3 = lracos(cost3);
        float t4 = lracos(dot(a3, b3));
        float cost7 = -dot(rcstackdir, b3);
        float t7 = lracos(cost7);
        float cost8 = dot(rcstackdir, a3);
        float t8 = lracos(cost8);

        float f2 = _f2(rcstackmod, CRST_F2, p);
        float f4t1 = _f4(t1, CRST_THETA1_T0, CRST_THETA1_TS, CRST_THETA1_TC, CRST_THETA1_A, CRST_THETA1_B);
        float f4t2 = _f4(t2, CRST_THETA2_T0, CRST_THETA2_TS, CRST_THETA2_TC, CRST_THETA2_A, CRST_THETA2_B);
        float f4t3 = _f4(t3, CRST_THETA3_T0, CRST_THETA3_TS, CRST_THETA3_TC, CRST_THETA3_A, CRST_THETA3_B);
        float f4t4 = _f4(t4, CRST_THETA4_T0, CRST_THETA4_TS, CRST_THETA4_TC, CRST_THETA4_A, CRST_THETA4_B)
                   + _f4(PI - t4, CRST_THETA4_T0, CRST_THETA4_TS, CRST_THETA4_TC, CRST_THETA4_A, CRST_THETA4_B);
        float f4t7 = _f4(t7, CRST_THETA7_T0, CRST_THETA7_TS, CRST_THETA7_TC, CRST_THETA7_A, CRST_THETA7_B)
                   + _f4(PI - t7, CRST_THETA7_T0, CRST_THETA7_TS, CRST_THETA7_TC, CRST_THETA7_A, CRST_THETA7_B);
        float f4t8 = _f4(t8, CRST_THETA8_T0, CRST_THETA8_TS, CRST_THETA8_TC, CRST_THETA8_A, CRST_THETA8_B)
                   + _f4(PI - t8, CRST_THETA8_T0, CRST_THETA8_TS, CRST_THETA8_TC, CRST_THETA8_A, CRST_THETA8_B);

        float cstk_energy = f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

        if(cstk_energy < 0.0f) {
            float f2D = _f2D(rcstackmod, CRST_F2, p);
            float f4t1D = -_f4D(t1, CRST_THETA1_T0, CRST_THETA1_TS, CRST_THETA1_TC, CRST_THETA1_A, CRST_THETA1_B);
            float f4t2D = -_f4D(t2, CRST_THETA2_T0, CRST_THETA2_TS, CRST_THETA2_TC, CRST_THETA2_A, CRST_THETA2_B);
            float f4t3D = _f4D(t3, CRST_THETA3_T0, CRST_THETA3_TS, CRST_THETA3_TC, CRST_THETA3_A, CRST_THETA3_B);
            float f4t4D = _f4D(t4, CRST_THETA4_T0, CRST_THETA4_TS, CRST_THETA4_TC, CRST_THETA4_A, CRST_THETA4_B)
                        - _f4D(PI - t4, CRST_THETA4_T0, CRST_THETA4_TS, CRST_THETA4_TC, CRST_THETA4_A, CRST_THETA4_B);
            float f4t7D = -_f4D(t7, CRST_THETA7_T0, CRST_THETA7_TS, CRST_THETA7_TC, CRST_THETA7_A, CRST_THETA7_B)
                        + _f4D(PI - t7, CRST_THETA7_T0, CRST_THETA7_TS, CRST_THETA7_TC, CRST_THETA7_A, CRST_THETA7_B);
            float f4t8D = _f4D(t8, CRST_THETA8_T0, CRST_THETA8_TS, CRST_THETA8_TC, CRST_THETA8_A, CRST_THETA8_B)
                        - _f4D(PI - t8, CRST_THETA8_T0, CRST_THETA8_TS, CRST_THETA8_TC, CRST_THETA8_A, CRST_THETA8_B);

            Ftmp = rcstackdir * (cstk_energy * f2D / f2);

            Ttmp -= stably_normalised(cross(a1, b1)) * (-cstk_energy * f4t1D / f4t1);

            Ftmp -= stably_normalised(b1 + rcstackdir * cost2) * (cstk_energy * f4t2D / (f4t2 * rcstackmod));

            float part = -cstk_energy * f4t3D / f4t3;
            Ftmp -= stably_normalised(a1 - rcstackdir * cost3) * (-part / rcstackmod);
            Ttmp += stably_normalised(cross(rcstackdir, a1)) * part;

            Ttmp -= stably_normalised(cross(a3, b3)) * (-cstk_energy * f4t4D / f4t4);

            Ftmp -= stably_normalised(b3 + rcstackdir * cost7) * (cstk_energy * f4t7D / (f4t7 * rcstackmod));

            part = -cstk_energy * f4t8D / f4t8;
            Ftmp -= stably_normalised(a3 - rcstackdir * cost8) * (-part / rcstackmod);
            Ttmp += stably_normalised(cross(rcstackdir, a3)) * part;

            Ttmp += cross(ppos_base, Ftmp);
            F += Ftmp;
        }
    }

    // ---- coaxial stacking ----
    float3 rstack = r + qpos_stack - ppos_stack;
    float rstackmodsqr = dot(rstack, rstack);
    if(SQR(CXST_RCLOW) < rstackmodsqr && rstackmodsqr < SQR(CXST_RCHIGH)) {
        float rstackmod = sqrt(rstackmodsqr);
        float3 rstackdir = rstack / rstackmod;

        float t1 = lracos(-dot(a1, b1));
        float t4 = lracos(dot(a3, b3));
        float cost5 = dot(a3, rstackdir);
        float t5 = lracos(cost5);
        float cost6 = -dot(b3, rstackdir);
        float t6 = lracos(cost6);

        float cosphi3 = 1.0f;
        float f5cosphi3 = 1.0f;
        float rbackrefmod = 1.0f;
        if(!use_ox2_cxst) {
            float3 rbackboneref = r + POS_BACK * b1 - POS_BACK * a1;
            rbackrefmod = length(rbackboneref);
            float3 rbackbonerefdir = rbackboneref / rbackrefmod;
            cosphi3 = dot(rstackdir, cross(rbackbonerefdir, a1));
            f5cosphi3 = _f5(cosphi3, CXST_F5_PHI3, p);
        }

        float f2 = _f2(rstackmod, CXST_F2, p);
        float f4t1 = use_ox2_cxst
            ? (_f4(t1, CXST_THETA1_T0_OXDNA2, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B)
               + _f4_pure_harmonic(t1, CXST_THETA1_SA, CXST_THETA1_SB))
            : (_f4(t1, CXST_THETA1_T0_OXDNA, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B)
               + _f4(2.0f * PI - t1, CXST_THETA1_T0_OXDNA, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B));
        float f4t4 = _f4(t4, CXST_THETA4_T0, CXST_THETA4_TS, CXST_THETA4_TC, CXST_THETA4_A, CXST_THETA4_B);
        float f4t5 = _f4(t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B)
                   + _f4(PI - t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B);
        float f4t6 = _f4(t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B)
                   + _f4(PI - t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B);

        float cxst_energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);

        if(cxst_energy < 0.0f) {
            float f2D = _f2D(rstackmod, CXST_F2, p);
            float f4t1D = use_ox2_cxst
                ? (-_f4D(t1, CXST_THETA1_T0_OXDNA2, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B)
                   - _f4D_pure_harmonic(t1, CXST_THETA1_SA, CXST_THETA1_SB))
                : (-_f4D(t1, CXST_THETA1_T0_OXDNA, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B)
                   + _f4D(2.0f * PI - t1, CXST_THETA1_T0_OXDNA, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B));
            float f4t4D = _f4D(t4, CXST_THETA4_T0, CXST_THETA4_TS, CXST_THETA4_TC, CXST_THETA4_A, CXST_THETA4_B);
            float f4t5D = _f4D(t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B)
                        - _f4D(PI - t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B);
            float f4t6D = -_f4D(t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B)
                        + _f4D(PI - t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B);

            Ftmp = rstackdir * (cxst_energy * f2D / f2);

            Ttmp -= stably_normalised(cross(a1, b1)) * (-cxst_energy * f4t1D / f4t1);
            Ttmp -= stably_normalised(cross(a3, b3)) * (-cxst_energy * f4t4D / f4t4);

            float part = cxst_energy * f4t5D / f4t5;
            Ftmp -= stably_normalised(a3 - rstackdir * cost5) / rstackmod * part;
            Ttmp -= stably_normalised(cross(rstackdir, a3)) * part;

            Ftmp -= stably_normalised(b3 + rstackdir * cost6) * (cxst_energy * f4t6D / (f4t6 * rstackmod));

            if(!use_ox2_cxst) {
                float f5cosphi3D = _f5D(cosphi3, CXST_F5_PHI3, p);
                float rbackrefmodcub = CUB(rbackrefmod);

                float a2b1 = dot(a2, b1);
                float a3b1 = dot(a3, b1);
                float ra1 = dot(rstackdir, a1);
                float ra2 = dot(rstackdir, a2);
                float ra3 = dot(rstackdir, a3);
                float rb1 = dot(rstackdir, b1);

                float parentesi = (ra3 * a2b1 - ra2 * a3b1);
                float dcdr    = -GAMMA * parentesi * (GAMMA * (ra1 - rb1) + rstackmod) / rbackrefmodcub;
                float dcda1b1 =  GAMMA * SQR(GAMMA) * parentesi / rbackrefmodcub;
                float dcda2b1 =  GAMMA * ra3 / rbackrefmod;
                float dcda3b1 = -GAMMA * ra2 / rbackrefmod;
                float dcdra1  = -SQR(GAMMA) * parentesi * rstackmod / rbackrefmodcub;
                float dcdra2  = -GAMMA * a3b1 / rbackrefmod;
                float dcdra3  =  GAMMA * a2b1 / rbackrefmod;
                float dcdrb1  = -dcdra1;

                part = cxst_energy * 2.0f * f5cosphi3D / f5cosphi3;

                Ftmp -= part * (rstackdir * dcdr
                                + ((a1 - rstackdir * ra1) * dcdra1
                                   + (a2 - rstackdir * ra2) * dcdra2
                                   + (a3 - rstackdir * ra3) * dcdra3
                                   + (b1 - rstackdir * rb1) * dcdrb1) / rstackmod);

                Ttmp += part * (cross(rstackdir, a1) * dcdra1
                                + cross(rstackdir, a2) * dcdra2
                                + cross(rstackdir, a3) * dcdra3);
                Ttmp -= part * (cross(a1, b1) * dcda1b1
                                + cross(a2, b1) * dcda2b1
                                + cross(a3, b1) * dcda3b1);
            }

            Ttmp += cross(ppos_stack, Ftmp);
            F += Ftmp;
        }
    }

    T += Ttmp;
}

// ---------------------------------------------------------------------------
//  Kernels
// ---------------------------------------------------------------------------
kernel void init_DNA_strand_ends(device int *is_strand_end [[buffer(0)]],
                                 device MetalBonds *bonds   [[buffer(1)]],
                                 constant InitStrandArgs &args [[buffer(2)]],
                                 uint2 tid [[thread_position_in_grid]]) {
    int idx = tid.x;
    if(idx >= args.N) return;
    MetalBonds b = bonds[idx];
    is_strand_end[idx] = (b.n3 == -1 || b.n5 == -1) ? 1 : 0;
}

kernel void dna_forces(device m_number4 *poss           [[buffer(0)]],
                       device m_number4 *orientations   [[buffer(1)]],
                       device m_number4 *forces         [[buffer(2)]],
                       device m_number4 *torques        [[buffer(3)]],
                       device int *matrix_neighs        [[buffer(4)]],
                       device int *number_neighs        [[buffer(5)]],
                       device MetalBonds *bonds         [[buffer(6)]],
                       constant DNAInteractionParams &params [[buffer(7)]],
                       constant MetalBox &box           [[buffer(8)]],
                       constant InitStrandArgs &args    [[buffer(9)]],
                       device float *energies           [[buffer(10)]],
                       uint2 tid [[thread_position_in_grid]]) {
    int idx = tid.x;
    if(idx >= args.N) return;

    if(energies) {
        for(int k = 0; k < 10; k++) energies[idx * 10 + k] = 0.0f;
    }

    float3 ppos = poss[idx].xyz;
    int ptype = (int) poss[idx].w;
    MetalBonds pb = bonds[idx];
    bool p_is_end = (pb.n3 == -1 || pb.n5 == -1);

    float3 a1, a2, a3;
    get_axes(orientations[idx], a1, a2, a3);

    float3 F = float3(0.0f);
    float3 T = float3(0.0f);

    // ---- bonded: 3' neighbour ----
    if(pb.n3 != -1) {
        int j = pb.n3;
        float3 qpos = poss[j].xyz;
        int qtype = (int) poss[j].w;
        float3 b1, b2, b3;
        get_axes(orientations[j], b1, b2, b3);

        float3 r = qpos - ppos;
        r = minimum_image(r, box);
        float3 dF = float3(0.0f);
        // qIsN3 = true : n5 partner = me (a*), n3 partner = j (b*)
        _bonded_part(true, r, ptype, a1, a2, a3, qtype, b1, b2, b3, dF, T, params);
        F += dF;
    }

    // ---- bonded: 5' neighbour ----
    if(pb.n5 != -1) {
        int j = pb.n5;
        float3 qpos = poss[j].xyz;
        int qtype = (int) poss[j].w;
        float3 b1, b2, b3;
        get_axes(orientations[j], b1, b2, b3);

        float3 r = ppos - qpos;
        r = minimum_image(r, box);
        float3 dF = float3(0.0f);
        // qIsN3 = false : n5 partner = j (b*), n3 partner = me (a*)
        _bonded_part(false, r, qtype, b1, b2, b3, ptype, a1, a2, a3, dF, T, params);
        F += dF;
    }

    // ---- non-bonded neighbours ----
    int n_neighs = number_neighs[idx];
    for(int i = 0; i < n_neighs; i++) {
        int j = matrix_neighs[i * args.N + idx];
        if(j == idx || j == pb.n3 || j == pb.n5) continue;

        float3 qpos = poss[j].xyz;
        int qtype = (int) poss[j].w;
        float3 r = qpos - ppos;
        r = minimum_image(r, box);

        // The Verlet list carries a skin margin; skip pairs that are beyond the
        // true interaction cutoff this step. (Matches the CPU pair filter.)
        if(dot(r, r) > params.sqr_rcut) continue;

        float3 b1, b2, b3;
        get_axes(orientations[j], b1, b2, b3);
        MetalBonds qb = bonds[j];
        bool q_is_end = (qb.n3 == -1 || qb.n5 == -1);

        float3 dF = float3(0.0f);
        _particle_particle_DNA_interaction(r, ptype, a1, a2, a3, qtype, b1, b2, b3,
                                           p_is_end, q_is_end, dF, T, params);
        F += dF;
    }

    // torque -> particle body frame (matches CUDA _vectors_transpose_c_number4_product)
    float3 Tbody = float3(dot(a1, T), dot(a2, T), dot(a3, T));

    forces[idx].xyz += F;
    torques[idx].xyz += Tbody;
}

// ---------------------------------------------------------------------------
//  Edge-list force kernels (Metal_list = edge)
//
//  Instead of one thread per particle recomputing every non-bonded pair twice
//  (once from each end), one thread per *edge* computes the pair once and
//  atomically distributes the force/torque to both partners. Torque is left in
//  the lab frame here; dna_forces_edge_bonded adds the bonded terms and does
//  the single lab -> body transform, exactly like the CUDA edge kernels.
// ---------------------------------------------------------------------------
kernel void dna_forces_edge_nonbonded(device m_number4 *poss              [[buffer(0)]],
                                      device m_number4 *orientations      [[buffer(1)]],
                                      device atomic_float *forces         [[buffer(2)]],
                                      device atomic_float *torques        [[buffer(3)]],
                                      device MetalEdgeBond *edge_list     [[buffer(4)]],
                                      device const int *n_edges           [[buffer(5)]],
                                      device MetalBonds *bonds            [[buffer(6)]],
                                      constant DNAInteractionParams &params [[buffer(7)]],
                                      constant MetalBox &box              [[buffer(8)]],
                                      uint2 tid [[thread_position_in_grid]]) {
    int e = tid.x;
    if(e >= n_edges[0]) return;

    MetalEdgeBond eb = edge_list[e];
    int pi = eb.from;
    int qi = eb.to;

    float3 ppos = poss[pi].xyz;
    float3 qpos = poss[qi].xyz;
    float3 r = qpos - ppos;
    r = minimum_image(r, box);
    if(dot(r, r) > params.sqr_rcut) return;

    int ptype = (int) poss[pi].w;
    int qtype = (int) poss[qi].w;

    float3 a1, a2, a3; get_axes(orientations[pi], a1, a2, a3);
    float3 b1, b2, b3; get_axes(orientations[qi], b1, b2, b3);

    MetalBonds pb = bonds[pi];
    MetalBonds qb = bonds[qi];
    bool p_is_end = (pb.n3 == -1 || pb.n5 == -1);
    bool q_is_end = (qb.n3 == -1 || qb.n5 == -1);

    float3 dF = float3(0.0f);
    float3 dT = float3(0.0f);
    _particle_particle_DNA_interaction(r, ptype, a1, a2, a3, qtype, b1, b2, b3,
                                       p_is_end, q_is_end, dF, dT, params);

    // particle p
    atomic_fetch_add_explicit(&forces[4 * pi + 0], dF.x, memory_order_relaxed);
    atomic_fetch_add_explicit(&forces[4 * pi + 1], dF.y, memory_order_relaxed);
    atomic_fetch_add_explicit(&forces[4 * pi + 2], dF.z, memory_order_relaxed);
    atomic_fetch_add_explicit(&torques[4 * pi + 0], dT.x, memory_order_relaxed);
    atomic_fetch_add_explicit(&torques[4 * pi + 1], dT.y, memory_order_relaxed);
    atomic_fetch_add_explicit(&torques[4 * pi + 2], dT.z, memory_order_relaxed);

    // particle q: Newton's 3rd law for the force, Allen Eq. 6 for the torque
    float3 tq = -dT + cross(r, dF);
    atomic_fetch_add_explicit(&forces[4 * qi + 0], -dF.x, memory_order_relaxed);
    atomic_fetch_add_explicit(&forces[4 * qi + 1], -dF.y, memory_order_relaxed);
    atomic_fetch_add_explicit(&forces[4 * qi + 2], -dF.z, memory_order_relaxed);
    atomic_fetch_add_explicit(&torques[4 * qi + 0], tq.x, memory_order_relaxed);
    atomic_fetch_add_explicit(&torques[4 * qi + 1], tq.y, memory_order_relaxed);
    atomic_fetch_add_explicit(&torques[4 * qi + 2], tq.z, memory_order_relaxed);
}

kernel void dna_forces_edge_bonded(device m_number4 *poss           [[buffer(0)]],
                                   device m_number4 *orientations   [[buffer(1)]],
                                   device m_number4 *forces         [[buffer(2)]],
                                   device m_number4 *torques        [[buffer(3)]],
                                   device MetalBonds *bonds         [[buffer(4)]],
                                   constant DNAInteractionParams &params [[buffer(5)]],
                                   constant MetalBox &box           [[buffer(6)]],
                                   constant InitStrandArgs &args    [[buffer(7)]],
                                   uint2 tid [[thread_position_in_grid]]) {
    int idx = tid.x;
    if(idx >= args.N) return;

    float3 ppos = poss[idx].xyz;
    int ptype = (int) poss[idx].w;
    MetalBonds pb = bonds[idx];

    float3 a1, a2, a3;
    get_axes(orientations[idx], a1, a2, a3);

    // start from the non-bonded accumulation left by dna_forces_edge_nonbonded
    float3 F = forces[idx].xyz;
    float3 T = torques[idx].xyz;

    if(pb.n3 != -1) {
        int j = pb.n3;
        int qtype = (int) poss[j].w;
        float3 b1, b2, b3; get_axes(orientations[j], b1, b2, b3);
        float3 r = poss[j].xyz - ppos;
        r = minimum_image(r, box);
        float3 dF = float3(0.0f);
        _bonded_part(true, r, ptype, a1, a2, a3, qtype, b1, b2, b3, dF, T, params);
        F += dF;
    }
    if(pb.n5 != -1) {
        int j = pb.n5;
        int qtype = (int) poss[j].w;
        float3 b1, b2, b3; get_axes(orientations[j], b1, b2, b3);
        float3 r = ppos - poss[j].xyz;
        r = minimum_image(r, box);
        float3 dF = float3(0.0f);
        _bonded_part(false, r, qtype, b1, b2, b3, ptype, a1, a2, a3, dF, T, params);
        F += dF;
    }

    float3 Tbody = float3(dot(a1, T), dot(a2, T), dot(a3, T));
    forces[idx].xyz = F;
    torques[idx].xyz = Tbody;
}
