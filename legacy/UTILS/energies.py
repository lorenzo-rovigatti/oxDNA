try:
    import numpy as np
except:
    import mynumpy as np
import base
from utils import get_angle

def _f1(r, r0, rlow, rhigh, rclow, rchigh, blow, bhigh, a, eps, shift):
    val = 0.
    if r < rchigh:
        if r > rhigh: val = eps * bhigh * (r - rchigh)**2
        elif r > rlow:
            tmp = 1. - exp(-(r - r0) * a)
            val = eps * tmp * tmp - shift
        elif r > rclow: val = eps * blow * (r - rclow)**2
            
    return val

def _f2(r, r0, rc, rlow, rhigh, rclow, rchigh, blow, bhigh, k):
    val = 0.
    if r < rchigh:
	    if r > rhigh: val = k * bhigh * (r - rchigh)**2
	    elif r > rlow: val = k * 0.5 * ((r - r0)**2 - (rc - r0)**2)
	    elif r > rclow: val = k * blow * (r - rclow)**2
    
    return val

def _f4(t, t0, ts, tc, a, b):
    val = 0.
    t -= t0
    if t < 0: t = -t
    
    if t < tc :
        if t > ts: val = b * (tc - t)**2
        else: val = 1. - a * t * t
        
    return val

def _f5(f, xs, xc, a, b):
    val = 0.
    
    if f > xc:
        if f < xs: val = b * (xc - f)**2
        elif f < 0.: val = 1. - a * f * f
        else: val = 1.

    return val

def _excluded_volume(rv, sigma, rstar, b, rc):
    energy = 0.
    r2 = np.dot(rv, rv)
    
    if r2 < rc*rc:
        if r2 > rstar*rstar:
            rrc = np.sqrt(r2) - rc
            energy = base.EXCL_EPS * b * rrc * rrc
        else:
            lj_part = sigma*sigma*sigma*sigma*sigma*sigma / (r2*r2*r2)
            energy = 4. * base.EXCL_EPS * (lj_part * lj_part - lj_part)
            
    return energy

def cross_stacking(p, q, box):
    rstack = q.pos_stack - p.pos_stack
    diff = box * np.rint (rstack / box)
    rstack -= diff
    rstackmod2 = np.dot(rstack, rstack)

    if (base.CRST_RCLOW*base.CRST_RCLOW) > rstackmod2 or rstackmod2 > (base.CRST_RCHIGH*base.CRST_RCHIGH): return 0

    rstackmod = np.sqrt(rstackmod2)
    rstackdir = rstack / rstackmod

    t1 = get_angle(-p._a1, q._a1)
    t2 = get_angle(-q._a1, rstackdir)
    t3 = get_angle(p._a1, rstackdir)
    t4 = get_angle(p._a3, q._a3)
    t7 = get_angle(-q._a3, rstackdir)
    t8 = get_angle(p._a3, rstackdir)

    f2 = _f2(rstackmod, base.CRST_R0, base.CRST_RC, base.CRST_RLOW, base.CRST_RHIGH, base.CRST_RCLOW, base.CRST_RCHIGH, base.CRST_BLOW, base.CRST_BHIGH, base.CRST_K)

    f4t1 = _f4(t1, base.CRST_THETA1_T0, base.CRST_THETA1_TS, base.CRST_THETA1_TC, base.CRST_THETA1_A, base.CRST_THETA1_B)
    f4t2 = _f4(t2, base.CRST_THETA2_T0, base.CRST_THETA2_TS, base.CRST_THETA2_TC, base.CRST_THETA2_A, base.CRST_THETA2_B)
    f4t3 = _f4(t3, base.CRST_THETA3_T0, base.CRST_THETA3_TS, base.CRST_THETA3_TC, base.CRST_THETA3_A, base.CRST_THETA3_B)
    f4t4 = _f4(t4, base.CRST_THETA4_T0, base.CRST_THETA4_TS, base.CRST_THETA4_TC, base.CRST_THETA4_A, base.CRST_THETA4_B) + _f4(np.pi - t4, base.CRST_THETA4_T0, base.CRST_THETA4_TS, base.CRST_THETA4_TC, base.CRST_THETA4_A, base.CRST_THETA4_B)
    f4t7 = _f4(t7, base.CRST_THETA7_T0, base.CRST_THETA7_TS, base.CRST_THETA7_TC, base.CRST_THETA7_A, base.CRST_THETA7_B) + _f4(np.pi - t7, base.CRST_THETA7_T0, base.CRST_THETA7_TS, base.CRST_THETA7_TC, base.CRST_THETA7_A, base.CRST_THETA7_B)
    f4t8 = _f4(t8, base.CRST_THETA8_T0, base.CRST_THETA8_TS, base.CRST_THETA8_TC, base.CRST_THETA8_A, base.CRST_THETA8_B) + _f4(np.pi - t8, base.CRST_THETA8_T0, base.CRST_THETA8_TS, base.CRST_THETA8_TC, base.CRST_THETA8_A, base.CRST_THETA8_B)
    
    return f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8

def coaxial(p, q, box):
    rstack = q.pos_stack - p.pos_stack
    diff = box * np.rint (rstack / box)
    rstack -= diff
    rstackmod2 = np.dot(rstack, rstack)

    if (base.CXST_RCLOW*base.CXST_RCLOW) > rstackmod2 or rstackmod2 > (base.CXST_RCHIGH*base.CXST_RCHIGH): return 0

    rstackmod = np.sqrt(rstackmod2)
    rstackdir = rstack / rstackmod

    t1 = get_angle(-p._a1, q._a1)
    t4 = get_angle( p._a3, q._a3)
    t5 = get_angle( p._a3, rstackdir)
    t6 = get_angle(-q._a3, rstackdir)

    rbackbone = q.pos_back - p.pos_back - diff
    cosphi3 = np.dot(rstackdir, (np.cross(rbackbone / np.sqrt(np.dot(rbackbone, rbackbone)), p._a1)))
    
    f2 = _f2(rstackmod, base.CXST_R0, base.CXST_RC, base.CXST_RLOW, base.CXST_RHIGH, base.CXST_RCLOW, base.CXST_RCHIGH, base.CXST_BLOW, base.CXST_BHIGH, base.CXST_K)
    f4t1 = _f4(t1, base.CXST_THETA1_T0, base.CXST_THETA1_TS, base.CXST_THETA1_TC, base.CXST_THETA1_A, base.CXST_THETA1_B) + _f4(2 * np.pi - t1, base.CXST_THETA1_T0, base.CXST_THETA1_TS, base.CXST_THETA1_TC, base.CXST_THETA1_A, base.CXST_THETA1_B)
    f4t4 = _f4(t4, base.CXST_THETA4_T0, base.CXST_THETA4_TS, base.CXST_THETA4_TC, base.CXST_THETA4_A, base.CXST_THETA4_B)
    f4t5 = _f4(t5, base.CXST_THETA5_T0, base.CXST_THETA5_TS, base.CXST_THETA5_TC, base.CXST_THETA5_A, base.CXST_THETA5_B) + _f4(np.pi - t5, base.CXST_THETA5_T0, base.CXST_THETA5_TS, base.CXST_THETA5_TC, base.CXST_THETA5_A, base.CXST_THETA5_B)
    f4t6 = _f4(t6, base.CXST_THETA6_T0, base.CXST_THETA6_TS, base.CXST_THETA6_TC, base.CXST_THETA6_A, base.CXST_THETA6_B) + _f4(np.pi - t6, base.CXST_THETA6_T0, base.CXST_THETA6_TS, base.CXST_THETA6_TC, base.CXST_THETA6_A, base.CXST_THETA6_B)
    f5cosphi3 = _f5(cosphi3, base.CXST_PHI3_XS, base.CXST_PHI3_XC, base.CXST_PHI3_A, base.CXST_PHI3_B)
    
    return f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi3

def excluded_volume(p, q, box, nearest=False):
    diff = box * np.rint((q.cm_pos - p.cm_pos) / box)

    # BASE-BASE
    rcenter = q.pos_base - p.pos_base - diff
    energy = _excluded_volume(rcenter, base.EXCL_S2, base.EXCL_R2, base.EXCL_B2, base.EXCL_RC2)
    
    # P-BASE vs. Q-BACK
    rcenter = q.pos_back - p.pos_base - diff
    energy += _excluded_volume(rcenter, base.EXCL_S3, base.EXCL_R3, base.EXCL_B3, base.EXCL_RC3)
    
    # P-BACK vs. Q-BASE
    rcenter = q.pos_base - p.pos_back - diff
    energy += _excluded_volume(rcenter, base.EXCL_S4, base.EXCL_R4, base.EXCL_B4, base.EXCL_RC4)
    
    if not nearest:
        # BACK-BACK
        rcenter = q.pos_back - p.pos_back - diff
        energy += _excluded_volume(rcenter, base.EXCL_S1, base.EXCL_R1, base.EXCL_B1, base.EXCL_RC1)

    return energy
