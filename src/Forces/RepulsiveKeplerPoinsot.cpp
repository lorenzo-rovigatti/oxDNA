#include "RepulsiveKeplerPoinsot.h"

#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

#include <cmath>
#include <vector>
#include <string>

static inline number _dot(const LR_vector &a, const LR_vector &b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

static inline LR_vector _cross(const LR_vector &a, const LR_vector &b) {
    return LR_vector(
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    );
}

static inline number _norm(const LR_vector &a) {
    return std::sqrt(_dot(a,a));
}

static inline LR_vector _normalize(const LR_vector &a) {
    const number n = _norm(a);
    if(n <= (number)0) return LR_vector(0,0,0);
    return LR_vector(a.x/n, a.y/n, a.z/n);
}

RepulsiveKeplerPoinsot::RepulsiveKeplerPoinsot() : BaseForce() {
    _centre      = LR_vector(0., 0., 0.);
    _rate        = (number)0.0;

    _apex        = (number)1.2;
    _base        = (number)0.7;
    _base_radius = (number)0.45;

    _kappa       = (number)25.0;
}

/**
 * Parses:
 *  - stiff (required)
 *  - rate  (optional, default 0)
 *  - apex (optional, default 1.2)
 *  - base (optional, default 0.7)
 *  - base_radius (optional, default 0.45)
 *  - kappa (optional, default 25)
 *  - center (optional, default 0,0,0)
 *  - particle (required): id or list or -1
 */
std::tuple<std::vector<int>, std::string> RepulsiveKeplerPoinsot::init(input_file &inp) {
    BaseForce::init(inp);

    getInputNumber(&inp, "stiff", &_stiff, 1);
    getInputNumber(&inp, "rate",  &_rate,  0);

    getInputNumber(&inp, "apex",        &_apex,        0);
    getInputNumber(&inp, "base",        &_base,        0);
    getInputNumber(&inp, "base_radius", &_base_radius, 0);
    getInputNumber(&inp, "kappa",       &_kappa,       0);

    if(_apex <= _base) throw oxDNAException("RepulsiveKeplerPoinsot: require apex > base");
    if(_base_radius <= (number)0.0) throw oxDNAException("RepulsiveKeplerPoinsot: require base_radius > 0");
    if(_kappa <= (number)0.0) _kappa = (number)25.0;

    std::string str;
    double tmpf[3];

    if(getInputString(&inp, "center", str, 0) == KEY_FOUND) {
        int tmpi = sscanf(str.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
        if(tmpi != 3) throw oxDNAException("Could not parse center %s", str.c_str());
        _centre = LR_vector(tmpf[0], tmpf[1], tmpf[2]);
    }

    std::string particles_string;
    getInputString(&inp, "particle", particles_string, 1);

    auto particle_ids =
        Utils::get_particles_from_string(
            CONFIG_INFO->particles(),
            particles_string,
            "RepulsiveKeplerPoinsot"
        );

    std::string desc = Utils::sformat(
        "RepulsiveKeplerPoinsot(stiff=%g, rate=%g, apex=%g, base=%g, base_radius=%g, kappa=%g, center=%g,%g,%g)",
        _stiff, _rate, _apex, _base, _base_radius, _kappa, _centre.x, _centre.y, _centre.z
    );

    return std::make_tuple(particle_ids, desc);
}

/**
 * Compute the "star metric" s(p) for a union of 12 pentagonal pyramidal spikes.
 *
 * For each spike with unit axis n:
 *   t = n·p
 *   q = p - t n  (perpendicular component)
 *   R(t) = base_radius * (apex - t)/(apex - base)  for t in [base, apex]
 *   support(q) = max over 5 directions in perpendicular plane (pentagon support)
 *   s = support(q) / R(t)
 *
 * Inside spike when: t in [base, apex] and s <= 1.
 * Union over spikes => pick smallest s among spikes (most "inside").
 *
 * We use a soft-max for support(q) to get a usable gradient for forces.
 */
static inline bool star_metric_and_grad(
    const LR_vector &p,          // already min-imaged displacement from centre
    number base, number apex, number base_radius,
    number kappa,
    number &s_best,
    LR_vector &grad_s_best
) {
    // 12 axes: icosahedron vertices (normalized)
    const number phi = (number)((1.0 + std::sqrt(5.0)) * 0.5);

    LR_vector axes[12] = {
        LR_vector(0,  1,  phi), LR_vector(0, -1,  phi),
        LR_vector(0,  1, -phi), LR_vector(0, -1, -phi),

        LR_vector( 1,  phi, 0), LR_vector(-1,  phi, 0),
        LR_vector( 1, -phi, 0), LR_vector(-1, -phi, 0),

        LR_vector( phi, 0,  1), LR_vector(-phi, 0,  1),
        LR_vector( phi, 0, -1), LR_vector(-phi, 0, -1)
    };

    // Pentagon directions in local (u,v) plane:
    const number two_pi = (number)(2.0 * M_PI);
    number cth[5], sth[5];
    for(int j=0;j<5;j++){
        const number th = two_pi * (number)j / (number)5.0;
        cth[j] = std::cos(th);
        sth[j] = std::sin(th);
    }

    bool found = false;
    s_best = (number)1e30;
    grad_s_best = LR_vector(0,0,0);

    for(int i=0;i<12;i++){
        const LR_vector n = _normalize(axes[i]);
        if(_norm(n) <= (number)0.0) continue;

        const number t = _dot(n, p);
        if(t < base || t > apex) continue;

        const number denom = (apex - base);
        if(denom <= (number)0.0) continue;

        const number R = base_radius * (apex - t) / denom;
        if(R <= (number)1e-12) continue;

        // Build orthonormal basis (u,v) perpendicular to n
        LR_vector a = (std::fabs(n.z) < (number)0.9) ? LR_vector(0,0,1) : LR_vector(0,1,0);
        LR_vector u = _normalize(_cross(n, a));
        if(_norm(u) <= (number)1e-12) continue;
        LR_vector v = _cross(n, u); // already unit if n,u unit and perpendicular

        // Perpendicular component q
        const LR_vector q = LR_vector(p.x - t*n.x, p.y - t*n.y, p.z - t*n.z);

        // Coordinates in (u,v)
        const number x = _dot(q, u);
        const number y = _dot(q, v);

        // Soft-max pentagon support:
        // dirs d_j = cosθ u + sinθ v
        // proj_j = x cosθ + y sinθ
        number proj[5];
        number m = -1e30;
        for(int j=0;j<5;j++){
            proj[j] = x*cth[j] + y*sth[j];
            if(proj[j] > m) m = proj[j];
        }

        // stable softmax weights
        number sumw = (number)0.0;
        number w[5];
        for(int j=0;j<5;j++){
            const number e = std::exp(kappa*(proj[j] - m));
            w[j] = e;
            sumw += e;
        }
        if(sumw <= (number)0.0) continue;

        // support = m + (1/kappa) log(sum exp(kappa*(proj-m)))
        const number support = m + (number)(1.0/kappa) * std::log(sumw);

        // gradient of support wrt (x,y): dS/dx = Σ softmax_j * cosθj, dS/dy = Σ softmax_j * sinθj
        number dSdx = (number)0.0, dSdy = (number)0.0;
        for(int j=0;j<5;j++){
            const number sj = w[j] / sumw;
            dSdx += sj * cth[j];
            dSdy += sj * sth[j];
        }

        // grad_p support = u*dSdx + v*dSdy
        const LR_vector grad_support = LR_vector(
            u.x*dSdx + v.x*dSdy,
            u.y*dSdx + v.y*dSdy,
            u.z*dSdx + v.z*dSdy
        );

        // s = support / R(t)
        const number s = support / R;

        if(s < s_best){
            // d(1/R)/dt
            // R(t) = base_radius*(apex - t)/denom => R'(t) = -base_radius/denom
            const number Rp = -base_radius / denom;
            const number dInvRdt = -(Rp) / (R*R); // d(1/R)/dt = -R'(t)/R^2

            // grad_p t = n, so grad_p(1/R) = d(1/R)/dt * n
            const LR_vector grad_invR = LR_vector(dInvRdt*n.x, dInvRdt*n.y, dInvRdt*n.z);

            // grad_p s = (1/R)*grad_support + support*grad_p(1/R)
            const number invR = (number)1.0 / R;
            const LR_vector grad_s = LR_vector(
                invR*grad_support.x + support*grad_invR.x,
                invR*grad_support.y + support*grad_invR.y,
                invR*grad_support.z + support*grad_invR.z
            );

            s_best = s;
            grad_s_best = grad_s;
            found = true;
        }
    }

    return found;
}

LR_vector RepulsiveKeplerPoinsot::value(llint step, LR_vector &pos) {
    const number growth = (number)1.0 + _rate * (number)step;
    if(growth <= (number)0.0) return LR_vector(0,0,0);

    const number base        = _base        * growth;
    const number apex        = _apex        * growth;
    const number base_radius = _base_radius * growth;

    if(apex <= base || base_radius <= (number)0.0) return LR_vector(0,0,0);

    // Min-image displacement from centre to particle
    LR_vector p = CONFIG_INFO->box->min_image(_centre, pos);

    number s;
    LR_vector grad_s;
    const bool ok = star_metric_and_grad(p, base, apex, base_radius, _kappa, s, grad_s);
    if(!ok) return LR_vector(0,0,0);

    // Apply repulsion only if inside the shape (s < 1)
    const number rc = (number)1.0;
    if(s >= rc) return LR_vector(0,0,0);

    // WCA-like in "s"
    const number x = (number)4;
    const number sigma = (number)1;
    const number epsilon = _stiff;

    const number s_safe = (s > (number)1e-9) ? s : (number)1e-9;

    const number A = std::pow(sigma / s_safe, x);

    // dU/ds = 4ε(2A - 1) dA/ds , dA/ds = -(x/s) A
    const number dUds = (number)4.0 * epsilon * ((number)2.0 * A - (number)1.0) * (-(x / s_safe) * A);

    // Force: F = - dU/ds * grad(s)
    const number scale = -dUds;

    return LR_vector(scale*grad_s.x, scale*grad_s.y, scale*grad_s.z);
}

number RepulsiveKeplerPoinsot::potential(llint step, LR_vector &pos) {
    const number growth = (number)1.0 + _rate * (number)step;
    if(growth <= (number)0.0) return (number)0.0;

    const number base        = _base        * growth;
    const number apex        = _apex        * growth;
    const number base_radius = _base_radius * growth;

    if(apex <= base || base_radius <= (number)0.0) return (number)0.0;

    LR_vector p = CONFIG_INFO->box->min_image(_centre, pos);

    number s;
    LR_vector grad_s_dummy;
    const bool ok = star_metric_and_grad(p, base, apex, base_radius, _kappa, s, grad_s_dummy);
    if(!ok) return (number)0.0;

    const number rc = (number)1.0;
    if(s >= rc) return (number)0.0;

    const number x = (number)4;
    const number sigma = (number)1;
    const number epsilon = _stiff;

    const number s_safe = (s > (number)1e-9) ? s : (number)1e-9;
    const number A = std::pow(sigma / s_safe, x);

    // U(s) = 4ε(A^2 - A) + ε  for s < 1
    const number U = (number)4.0 * epsilon * (A*A - A) + epsilon;
    return U;
}
