#include "StressInBox.h"

#include <cmath>
#include <sstream>
#include <iomanip>
#include <fstream>

// Important: use relative includes from src/Observables
#include "../Boxes/BaseBox.h"
#include "../Particles/BaseParticle.h"
// For subset stress helper
#include "../Interactions/BaseInteraction.h"

StressInBox::StressInBox() : BaseObservable() {}

StressInBox::~StressInBox() {
    if(_dbg.is_open()) _dbg.close();
}

bool StressInBox::require_data_on_CPU() {
    return true;
}

void StressInBox::get_settings(input_file &my_inp, input_file &sim_inp) {
    BaseObservable::get_settings(my_inp, sim_inp);

    // Side length (optional)
    getInputDouble(&my_inp, "L", &_L, 0);

    // Center (optional; default 0,0,0)
    double cx = 0.0, cy = 0.0, cz = 0.0;
    getInputDouble(&my_inp, "center_x", &cx, 0);
    getInputDouble(&my_inp, "center_y", &cy, 0);
    getInputDouble(&my_inp, "center_z", &cz, 0);
    _center = LR_vector(cx, cy, cz);

    // Debug controls (all optional)
    getInputBool(&my_inp, "debug", &_debug, 0);
    getInputLLInt(&my_inp, "debug_every", &_debug_every, 0);
    getInputInt(&my_inp, "debug_samples", &_debug_samples, 0);
    getInputString(&my_inp, "debug_file", _debug_filename, 0);

    // Precompute
    _halfL = 0.5 * _L;
    _Vbox  = _L * _L * _L;
    if(_Vbox <= 0.0) {
        // Avoid oxDNAException include here; just clamp to prevent div0
        _Vbox = 1.0;
    }
}

void StressInBox::init() {
    BaseObservable::init();

    if (_debug) {
        bool empty = true;
        std::ifstream fin(_debug_filename.c_str());
        if (fin.good())
            empty = (fin.peek() == std::ifstream::traits_type::eof());

        _dbg.open(_debug_filename.c_str(), std::ios::out | std::ios::app);
        _dbg.setf(std::ios::scientific);
        _dbg << std::setprecision(16);

        if (empty) {
            _dbg << "step n_in "
                 << "Fsum_abs Fmax_abs nFpos "
                 << "dr_abs_sum dr_abs_max "
                 << "virial_abs_sum "
                 << "xx xy xz yx yy yz zx zy zz\n";
            _dbg.flush();
        }
    }
}

void StressInBox::update_data(llint curr_step) {
    // reset tensor
    xx = yy = zz = 0.0;
    xy = xz = yz = 0.0;
    yx = zx = zy = 0.0;

    // reset diagnostics
    n_in = 0;
    Fsum_abs = 0.0;
    Fmax_abs = 0.0;
    nFpos = 0;

    Fsum_vec = LR_vector(0.,0.,0.);
    dr_abs_sum = 0.0;
    dr_abs_max = 0.0;
    virial_abs_sum = 0.0;

    BaseBox *box = _config_info->box;

    // NOTE: in your fork, particles is a METHOD, not a member
    auto &parts = _config_info->particles();
    const int N = (int) parts.size();

    int dumped = 0;

    // Build the subset once and compute the stress using BaseInteraction.
    // We still compute diagnostics in the same loop to help identify cancellation.
    std::vector<BaseParticle *> subset;
    subset.reserve(N);

    for(int i = 0; i < N; i++) {
        BaseParticle *p = parts[i];

        // displacement from center with minimum image (PBC aware)
        LR_vector dr = box->min_image(p->pos, _center);

        if(!_inside_box(dr)) continue;

        n_in++;

        // total force on particle
        LR_vector F = p->force;

        const double fmod = F.module();
        Fsum_abs += fmod;
        if(fmod > Fmax_abs) Fmax_abs = fmod;
        if(fmod > 0.0) nFpos++;

        Fsum_vec += F;

        const double dr_abs = dr.module();
        dr_abs_sum += dr_abs;
        if(dr_abs > dr_abs_max) dr_abs_max = dr_abs;

        subset.push_back(p);

        // cancellation-proof magnitude
        virial_abs_sum += std::fabs(dr.x*F.x) + std::fabs(dr.x*F.y) + std::fabs(dr.x*F.z)
                       +  std::fabs(dr.y*F.x) + std::fabs(dr.y*F.y) + std::fabs(dr.y*F.z)
                       +  std::fabs(dr.z*F.x) + std::fabs(dr.z*F.y) + std::fabs(dr.z*F.z);

        // --- debug sample ---
        if (_debug &&
            _debug_every > 0 &&
            (curr_step % _debug_every == 0) &&
            dumped < _debug_samples)
        {
            _dbg << "# sample step=" << curr_step
                 << " i=" << i
                 << " dr=(" << dr.x << "," << dr.y << "," << dr.z << ")"
                 << " F=("  << F.x  << "," << F.y  << "," << F.z  << ")"
                 << " |F|=" << fmod << "\n";
            dumped++;
        }
    }

    // Compute the tensor.
    // If the box contains ALL particles, match StressAutocorrelation exactly by using
    // compute_standard_stress_tensor()/stress_tensor(). Otherwise, fall back to the
    // subset estimator (virial from total forces, normalised by Vbox).
    if(_config_info->interaction != NULL) {
        const int N_total = (int) parts.size();
        if(n_in == N_total) {
            if(!_config_info->interaction->has_custom_stress_tensor()) {
                _config_info->interaction->compute_standard_stress_tensor();
            }
            StressTensor st = _config_info->interaction->stress_tensor();
            xx = st[0];
            yy = st[1];
            zz = st[2];
            xy = st[3];
            xz = st[4];
            yz = st[5];
        } else {
            StressTensor st = _config_info->interaction->stress_tensor_subset(
                subset, _center, (number) _Vbox, true /* use_min_image */, false /* include_kinetic */
            );
            xx = st[0];
            yy = st[1];
            zz = st[2];
            xy = st[3];
            xz = st[4];
            yz = st[5];
        }
        yx = xy;
        zx = xz;
        zy = yz;
    }

    if(_debug && _dbg.is_open()) {
        _dbg << curr_step << " " << n_in << " "
             << Fsum_abs << " "
             << Fsum_vec.x << " " << Fsum_vec.y << " " << Fsum_vec.z << " "
             << dr_abs_sum << " " << dr_abs_max << " "
             << virial_abs_sum << " "
             << xx << " " << yy << " " << zz << " "
             << xy << " " << xz << " " << yz
             << "\n";
        _dbg.flush();
    }
}

std::string StressInBox::get_output_string(llint curr_step) {
    std::ostringstream out;
    out.setf(std::ios::scientific);
    out << std::setprecision(4);

    // Main observable output line
    out << curr_step << " "
        << xx << " " << yy << " " << zz << " "
        << xy << " " << xz << " " << yz << " "
        << n_in << " "
        << Fsum_abs << " " << Fmax_abs << " " << nFpos << " "
        << Fsum_vec.x << " " << Fsum_vec.y << " " << Fsum_vec.z << " "
        << dr_abs_sum << " " << dr_abs_max << " "
        << virial_abs_sum;
    return out.str();
}