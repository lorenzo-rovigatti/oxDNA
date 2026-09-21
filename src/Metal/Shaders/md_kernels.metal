/**
 * @file    md_kernels.metal
 * @brief   Metal kernels for molecular dynamics
 *
 * Core MD kernels equivalent to CUDA MD kernels
 */

#include <metal_stdlib>
#include "common.metal"
#include "df64.h"
using namespace metal;

/**
 * @brief Update particle velocities (first half of velocity Verlet)
 */
inline m_number4 quat_multiply_md(m_number4 a, m_number4 b) {
    return m_number4(
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z
    );
}

inline m_number4 update_orientation_from_L(m_number4 L, m_number4 old_orientation, m_number dt) {
    m_number3 Lvec = L.xyz;
    m_number norm = length(Lvec);
    if(norm <= (m_number) 1e-12f) {
        return old_orientation;
    }

    m_number3 axis = Lvec / norm;
    m_number half_theta = dt * norm * (m_number) 0.5f;
    m_number s = sin(half_theta);
    m_number c = cos(half_theta);

    m_number4 R = m_number4(axis.x * s, axis.y * s, axis.z * s, c);
    m_number4 updated = quat_multiply_md(old_orientation, R);

    m_number qnorm = length(updated);
    if(qnorm <= (m_number) 1e-12f) {
        return old_orientation;
    }
    return updated / qnorm;
}

kernel void first_step_velocity_verlet(
    device m_number4 *positions [[buffer(0)]],
    device m_number4 *orientations [[buffer(1)]],
    device m_number4 *velocities [[buffer(2)]],
    device m_number4 *angular_momenta [[buffer(3)]],
    device m_number4 *forces [[buffer(4)]],
    device m_number4 *torques [[buffer(5)]],
    constant MetalBox &box [[buffer(6)]],
    constant m_number &dt [[buffer(7)]],
    constant m_number &dt_half [[buffer(8)]],
    constant int &N [[buffer(9)]],
    uint gid [[thread_position_in_grid]])
{
    if((int) gid >= N) {
        return;
    }

    // v += F * dt / (2m)
    m_number4 v = velocities[gid];
    const m_number4 F = forces[gid];
    v.x += F.x * dt_half;
    v.y += F.y * dt_half;
    v.z += F.z * dt_half;

    // r += v * dt
    m_number4 r = positions[gid];
    r.x += v.x * dt;
    r.y += v.y * dt;
    r.z += v.z * dt;

    m_number4 L = angular_momenta[gid];
    const m_number4 T = torques[gid];
    L.x += T.x * dt_half;
    L.y += T.y * dt_half;
    L.z += T.z * dt_half;

    positions[gid] = r;
    velocities[gid] = v;
    angular_momenta[gid] = L;
    orientations[gid] = update_orientation_from_L(L, orientations[gid], dt);
}

/**
 * @brief Update velocities (second half of velocity Verlet)
 */
kernel void second_step_velocity_verlet(
    device m_number4 *velocities [[buffer(0)]],
    device m_number4 *angular_momenta [[buffer(1)]],
    device m_number4 *forces [[buffer(2)]],
    device m_number4 *torques [[buffer(3)]],
    constant m_number &dt_half [[buffer(4)]],
    constant int &N [[buffer(5)]],
    uint gid [[thread_position_in_grid]])
{
    if((int) gid >= N) {
        return;
    }

    // v += F * dt / (2m)
    m_number4 v = velocities[gid];
    m_number4 F = forces[gid];
    v.x += F.x * dt_half;
    v.y += F.y * dt_half;
    v.z += F.z * dt_half;
    v.w = (v.x * v.x + v.y * v.y + v.z * v.z) * (m_number) 0.5f;
    velocities[gid] = v;

    m_number4 L = angular_momenta[gid];
    m_number4 T = torques[gid];
    L.x += T.x * dt_half;
    L.y += T.y * dt_half;
    L.z += T.z * dt_half;
    L.w = (L.x * L.x + L.y * L.y + L.z * L.z) * (m_number) 0.5f;
    angular_momenta[gid] = L;
}

/**
 * @brief Update orientations using quaternions
 */
kernel void update_orientations(
    device m_number4 *orientations [[buffer(0)]],
    device m_number4 *angular_velocities [[buffer(1)]],
    device m_number4 *torques [[buffer(2)]],
    constant m_number &dt [[buffer(3)]],
    uint gid [[thread_position_in_grid]])
{
    // Get current orientation quaternion
    m_number4 q = orientations[gid];

    // Get angular velocity
    m_number3 omega = m_number3(angular_velocities[gid].x,
                                angular_velocities[gid].y,
                                angular_velocities[gid].z);

    // Quaternion derivative: dq/dt = 0.5 * q * omega_quat (Body frame)
    m_number4 omega_quat = m_number4(omega.x, omega.y, omega.z, 0.0);

    // Quaternion multiplication q * omega
    // Real: -v.w
    // Vec: s*w + v x w
    // q = (x,y,z,w). w=scalar. Note: my q struct has w at end?
    // q.w is scalar part (based get_axes usage).
    // q.x, q.y, q.z is vector part.
    // omega = (Lx, Ly, Lz).
    // dq.x = q.w * Lx + (q.y * Lz - q.z * Ly)
    // dq.y = q.w * Ly + (q.z * Lx - q.x * Lz)
    // dq.z = q.w * Lz + (q.x * Ly - q.y * Lx)
    // dq.w = - (q.x * Lx + q.y * Ly + q.z * Lz)
    
    m_number4 dq;
    dq.x = 0.5 * (q.w * omega.x + q.y * omega.z - q.z * omega.y);
    dq.y = 0.5 * (q.w * omega.y + q.z * omega.x - q.x * omega.z);
    dq.z = 0.5 * (q.w * omega.z + q.x * omega.y - q.y * omega.x);
    dq.w = -0.5 * (q.x * omega.x + q.y * omega.y + q.z * omega.z);

    // Update quaternion
    q += dq * dt;

    // Normalize quaternion
    m_number norm = sqrt(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
    orientations[gid] = q / norm;
}

/**
 * @brief Update angular momenta (L += T * dt)
 */
kernel void update_angular_momenta(
    device m_number4 *angular_momenta [[buffer(0)]],
    device m_number4 *torques [[buffer(1)]],
    constant m_number &dt [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    // L += T * dt
    angular_momenta[gid].x += torques[gid].x * dt;
    angular_momenta[gid].y += torques[gid].y * dt;
    angular_momenta[gid].z += torques[gid].z * dt;
}

/**
 * @brief Compute kinetic energy
 */
kernel void compute_kinetic_energy(
    device m_number4 *velocities [[buffer(0)]],
    device m_number4 *positions [[buffer(1)]],
    device m_number *kinetic_energy [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    m_number mass = positions[gid].w;
    m_number3 v = m_number3(velocities[gid].x, velocities[gid].y, velocities[gid].z);
    m_number v_sqr = dot(v, v);

    kinetic_energy[gid] = 0.5 * mass * v_sqr;
}

/**
 * @brief Zero forces array
 */
kernel void zero_forces(
    device m_number4 *forces [[buffer(0)]],
    uint gid [[thread_position_in_grid]])
{
    forces[gid] = m_number4(0.0, 0.0, 0.0, 0.0);
}

/**
 * @brief Zero torques array
 */
kernel void zero_torques(
    device m_number4 *torques [[buffer(0)]],
    uint gid [[thread_position_in_grid]])
{
    torques[gid] = m_number4(0.0, 0.0, 0.0, 0.0);
}

/**
 * @brief Copy buffer
 */
kernel void copy_buffer_m_number4(
    device m_number4 *dest [[buffer(0)]],
    device m_number4 *src [[buffer(1)]],
    uint gid [[thread_position_in_grid]])
{
    dest[gid] = src[gid];
}

/**
 * @brief Reduction sum for m_number4
 */
kernel void reduce_sum_m_number4(
    device m_number4 *input [[buffer(0)]],
    device m_number4 *output [[buffer(1)]],
    threadgroup m_number4 *shared [[threadgroup(0)]],
    uint gid [[thread_position_in_grid]],
    uint tid [[thread_position_in_threadgroup]],
    uint threads [[threads_per_threadgroup]])
{
    // Load input into shared memory
    shared[tid] = input[gid];
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // Parallel reduction
    for(uint s = threads / 2; s > 0; s >>= 1) {
        if(tid < s) {
            shared[tid] += shared[tid + s];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }

    // Write result
    if(tid == 0) {
        output[0] = shared[0];
    }
}

// ============================================================================
// Mixed-precision (double-float, df64) velocity-Verlet integration.
//
// `mixed` keeps positions, velocities and angular momenta as df64 pairs
// (hi/lo float4 buffers) while forces/torques stay float32 — the Metal
// analogue of the CUDA mixed backend, which uses hardware double for the
// same quantities. The float mirror buffers (positions/velocities/momenta)
// are kept in sync every step: force kernels, Verlet-list checks and
// thermostats all read the mirrors. Orientations stay float32 quaternions
// (renormalized every step, exactly as in the float path).
// ============================================================================

struct df64_4 {
    float4 hi;
    float4 lo;
};

inline df64_4 df4_from_f4(float4 a) {
    df64_4 r;
    r.hi = a;
    r.lo = float4(0.0f);
    return r;
}

inline float4 df4_to_f4(df64_4 a) {
    return a.hi + a.lo;
}

inline df64_4 df4_add_df4(df64_4 a, df64_4 b) {
    df64_4 r;
    df64 cx = df_add_df(df64(a.hi.x, a.lo.x), df64(b.hi.x, b.lo.x));
    df64 cy = df_add_df(df64(a.hi.y, a.lo.y), df64(b.hi.y, b.lo.y));
    df64 cz = df_add_df(df64(a.hi.z, a.lo.z), df64(b.hi.z, b.lo.z));
    r.hi = float4(cx.hi, cy.hi, cz.hi, a.hi.w);
    r.lo = float4(cx.lo, cy.lo, cz.lo, 0.0f);
    return r;
}

inline df64_4 df4_mul_f(df64_4 a, float b) {
    df64_4 r;
    df64 cx = df_mul_f(df64(a.hi.x, a.lo.x), b);
    df64 cy = df_mul_f(df64(a.hi.y, a.lo.y), b);
    df64 cz = df_mul_f(df64(a.hi.z, a.lo.z), b);
    r.hi = float4(cx.hi, cy.hi, cz.hi, a.hi.w);
    r.lo = float4(cx.lo, cy.lo, cz.lo, 0.0f);
    return r;
}

/**
 * @brief First half of velocity Verlet in df64: v += F*dt/2, r += v*dt.
 */
kernel void first_step_mixed(
    device float4 *poss [[buffer(0)]],
    device float4 *orientations [[buffer(1)]],
    device float4 *vels_mir [[buffer(2)]],
    device float4 *ls_mir [[buffer(3)]],
    device float4 *forces [[buffer(4)]],
    device float4 *torques [[buffer(5)]],
    device float4 *poss_hi [[buffer(6)]],
    device float4 *poss_lo [[buffer(7)]],
    device float4 *vels_hi [[buffer(8)]],
    device float4 *vels_lo [[buffer(9)]],
    device float4 *ls_hi [[buffer(10)]],
    device float4 *ls_lo [[buffer(11)]],
    constant float &dt [[buffer(12)]],
    constant float &dt_half [[buffer(13)]],
    constant int &N [[buffer(14)]],
    uint gid [[thread_position_in_grid]])
{
    if((int) gid >= N) {
        return;
    }

    // v += F * dt/2 in double-float (exact product via TwoProd+FMA)
    df64_4 v;
    v.hi = vels_hi[gid];
    v.lo = vels_lo[gid];
    df64_4 dv = df4_mul_f(df4_from_f4(forces[gid]), dt_half);
    v = df4_add_df4(v, dv);
    vels_hi[gid] = v.hi;
    vels_lo[gid] = v.lo;
    float4 vf = df4_to_f4(v);
    vels_mir[gid] = float4(vf.x, vf.y, vf.z, 0.0f);

    // r += v * dt in double-float
    df64_4 r;
    r.hi = poss_hi[gid];
    r.lo = poss_lo[gid];
    r = df4_add_df4(r, df4_mul_f(v, dt));
    poss_hi[gid] = r.hi;
    poss_lo[gid] = r.lo;
    float4 rf = df4_to_f4(r);
    poss[gid] = float4(rf.x, rf.y, rf.z, poss[gid].w);

    // L += T * dt/2 in double-float
    df64_4 L;
    L.hi = ls_hi[gid];
    L.lo = ls_lo[gid];
    df64_4 dL = df4_mul_f(df4_from_f4(torques[gid]), dt_half);
    L = df4_add_df4(L, dL);
    ls_hi[gid] = L.hi;
    ls_lo[gid] = L.lo;
    float4 Lf = df4_to_f4(L);
    ls_mir[gid] = float4(Lf.x, Lf.y, Lf.z, 0.0f);

    orientations[gid] = update_orientation_from_L(Lf, orientations[gid], dt);
}

/**
 * @brief Second half of velocity Verlet in df64: v += F*dt/2.
 */
kernel void second_step_mixed(
    device float4 *vels_mir [[buffer(0)]],
    device float4 *ls_mir [[buffer(1)]],
    device float4 *forces [[buffer(2)]],
    device float4 *torques [[buffer(3)]],
    device float4 *vels_hi [[buffer(4)]],
    device float4 *vels_lo [[buffer(5)]],
    device float4 *ls_hi [[buffer(6)]],
    device float4 *ls_lo [[buffer(7)]],
    constant float &dt_half [[buffer(8)]],
    constant int &N [[buffer(9)]],
    uint gid [[thread_position_in_grid]])
{
    if((int) gid >= N) {
        return;
    }

    df64_4 v;
    v.hi = vels_hi[gid];
    v.lo = vels_lo[gid];
    v = df4_add_df4(v, df4_mul_f(df4_from_f4(forces[gid]), dt_half));
    vels_hi[gid] = v.hi;
    vels_lo[gid] = v.lo;
    float4 vf = df4_to_f4(v);
    vels_mir[gid] = float4(vf.x, vf.y, vf.z,
                           (vf.x * vf.x + vf.y * vf.y + vf.z * vf.z) * 0.5f);

    df64_4 L;
    L.hi = ls_hi[gid];
    L.lo = ls_lo[gid];
    L = df4_add_df4(L, df4_mul_f(df4_from_f4(torques[gid]), dt_half));
    ls_hi[gid] = L.hi;
    ls_lo[gid] = L.lo;
    float4 Lf = df4_to_f4(L);
    ls_mir[gid] = float4(Lf.x, Lf.y, Lf.z,
                         (Lf.x * Lf.x + Lf.y * Lf.y + Lf.z * Lf.z) * 0.5f);
}

/**
 * @brief Re-split the float velocity mirrors into df64 after the GPU
 * thermostat has rescaled them (hi = mirror, lo = 0).
 */
kernel void mixed_sync_df_vels(
    device float4 *vels_mir [[buffer(0)]],
    device float4 *ls_mir [[buffer(1)]],
    device float4 *vels_hi [[buffer(2)]],
    device float4 *vels_lo [[buffer(3)]],
    device float4 *ls_hi [[buffer(4)]],
    device float4 *ls_lo [[buffer(5)]],
    constant int &N [[buffer(6)]],
    uint gid [[thread_position_in_grid]])
{
    if((int) gid >= N) {
        return;
    }

    float4 v = vels_mir[gid];
    vels_hi[gid] = float4(v.x, v.y, v.z, 0.0f);
    vels_lo[gid] = float4(0.0f);
    float4 L = ls_mir[gid];
    ls_hi[gid] = float4(L.x, L.y, L.z, 0.0f);
    ls_lo[gid] = float4(0.0f);
}
