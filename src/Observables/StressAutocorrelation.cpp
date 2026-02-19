// StressAutocorrelation.cpp

#include "StressAutocorrelation.h"

#include <sstream>
#include <fstream>   // ADD
#include <sys/stat.h> // optional; if unavailable, remove and use simple ifstream test

StressAutocorrelation::StressAutocorrelation() :
                BaseObservable() {

}

StressAutocorrelation::~StressAutocorrelation() {
    if(_tensor_out.is_open()) {
        _tensor_out.close();
    }
    if(_tensor_cyl_out.is_open()) _tensor_cyl_out.close();
}

void StressAutocorrelation::serialise() {
    _sigma_P->serialise("sigma_P.dat");
    _sigma_xx->serialise("sigma_xx.dat");
    _sigma_yy->serialise("sigma_yy.dat");
    _sigma_zz->serialise("sigma_zz.dat");
    _sigma_xy->serialise("sigma_xy.dat");
    _sigma_yz->serialise("sigma_yz.dat");
    _sigma_xz->serialise("sigma_xz.dat");

    _N_xy->serialise("N_xy.dat");
    _N_yz->serialise("N_yz.dat");   // FIX
    _N_xz->serialise("N_xz.dat");   // FIX (was N_zx.dat)
}

std::shared_ptr<StressAutocorrelation::Level> StressAutocorrelation::_deserialise(std::string filename) {
    std::shared_ptr<Level> res = std::make_shared<Level>(_m, _p, 0);

    std::ifstream inp(filename);
    if(inp) {
        llint step;

        inp.ignore(32768, '=');
        inp >> step;
        if(step != CONFIG_INFO->curr_step) {
            inp.close();
            throw oxDNAException("StressAutocorrelation: the timestep found in the '%s' file (%lld) "
                    "does not match the one read by the initial configuration (%lld). Check that "
                    "the initial configuration file is correct and that restart_step_counter = false",
                    filename.c_str(), step, CONFIG_INFO->curr_step);
        }

        res->load_from_file(inp);
        inp.close();
    }
    else {
        OX_LOG(Logger::LOG_WARNING, "StressAutocorrelation: file '%s' not found", filename.c_str());
    }

    return res;
}

void StressAutocorrelation::get_settings(input_file &my_inp, input_file &sim_inp) {
    BaseObservable::get_settings(my_inp, sim_inp);

    if(_update_every == 0) {
        throw oxDNAException("StressAutocorrelation: update_every should be larger than 0");
    }

    getInputUInt(&my_inp, "m", &_m, 0);
    getInputUInt(&my_inp, "p", &_p, 0);
    getInputBool(&my_inp, "anisotropic", &_anisotropic, 0);

    getInputDouble(&sim_inp, "dt", &_delta_t, 1);
    _delta_t *= _update_every;

    getInputBool(&my_inp, "serialise", &_enable_serialisation, 0);

    // =========================
    // NEW optional inputs
    // =========================
    getInputBool(&my_inp, "dump_tensor", &_dump_tensor, 0);
    getInputString(&my_inp, "tensor_file", _tensor_filename, 0);

    getInputBool(&my_inp, "dump_tensor_cyl", &_dump_tensor_cyl, 0);
    getInputString(&my_inp, "tensor_cyl_file", _tensor_cyl_filename, 0);

    // Optional origin (if not provided, default to box center in x/y during init)
    double ox, oy, oz;
    if(getInputDouble(&my_inp, "cyl_origin_x", &ox, 0) == KEY_FOUND &&
    getInputDouble(&my_inp, "cyl_origin_y", &oy, 0) == KEY_FOUND) {
        getInputDouble(&my_inp, "cyl_origin_z", &oz, 0); // optional
        _cyl_origin = LR_vector(ox, oy, oz);
        _cyl_origin_set = true;
    }

    OX_LOG(Logger::LOG_INFO,
        "Initialising a StressAutocorrelation obs with m = %d, p = %d, anisotropic = %d, dump_tensor = %d (%s)",
        _m, _p, _anisotropic, _dump_tensor, _tensor_filename.c_str()
    );
}

static bool _file_has_content(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if(!f) return false;
    f.seekg(0, std::ios::end);
    return (f.tellg() > 0);
}

void StressAutocorrelation::init() {
    if(_enable_serialisation) {
        _sigma_P  = _deserialise("sigma_P.dat");
        _sigma_xx = _deserialise("sigma_xx.dat");
        _sigma_yy = _deserialise("sigma_yy.dat");
        _sigma_zz = _deserialise("sigma_zz.dat");
        _sigma_xy = _deserialise("sigma_xy.dat");
        _sigma_yz = _deserialise("sigma_yz.dat");
        _sigma_xz = _deserialise("sigma_xz.dat");

        _N_xy = _deserialise("N_xy.dat");
        _N_yz = _deserialise("N_yz.dat");
        _N_xz = _deserialise("N_xz.dat"); // FIX (was N_zx.dat)
    }
    else {
        _sigma_P  = std::make_shared<Level>(_m, _p, 0);
        _sigma_xx = std::make_shared<Level>(_m, _p, 0);
        _sigma_yy = std::make_shared<Level>(_m, _p, 0);
        _sigma_zz = std::make_shared<Level>(_m, _p, 0);
        _sigma_xy = std::make_shared<Level>(_m, _p, 0);
        _sigma_yz = std::make_shared<Level>(_m, _p, 0);
        _sigma_xz = std::make_shared<Level>(_m, _p, 0);
        _N_xy     = std::make_shared<Level>(_m, _p, 0);
        _N_yz     = std::make_shared<Level>(_m, _p, 0);
        _N_xz     = std::make_shared<Level>(_m, _p, 0);
    }

    // =========================
    // NEW: open tensor file in APPEND mode
    // =========================
    if(_dump_tensor) {
        const bool write_header = !_file_has_content(_tensor_filename);

        _tensor_out.open(_tensor_filename, std::ios::out | std::ios::app);
        if(!_tensor_out.good()) {
            throw oxDNAException("StressAutocorrelation: cannot open tensor_file '%s' for appending",
                                 _tensor_filename.c_str());
        }

        if(write_header) {
            _tensor_out
                << "step  "
                << "xx xy xz  "
                << "yx yy yz  "
                << "zx zy zz"
                << std::endl;
        }
    }
    if(_dump_tensor_cyl) {
        // Default origin: box center in x/y if not provided
        if(!_cyl_origin_set) {
            LR_vector sides = _config_info->box->box_sides();
            _cyl_origin = LR_vector(0.5*sides.x, 0.5*sides.y, 0.0);
        }

        const bool write_header = !_file_has_content(_tensor_cyl_filename);
        _tensor_cyl_out.open(_tensor_cyl_filename, std::ios::out | std::ios::app);
        if(!_tensor_cyl_out.good()) {
            throw oxDNAException("StressAutocorrelation: cannot open tensor_cyl_file '%s' for appending",
                                _tensor_cyl_filename.c_str());
        }

        if(write_header) {
            _tensor_cyl_out
                << "# step  "
                << "rr rt rz  "
                << "tr tt tz  "
                << "zr zt zz"
                << std::endl;
        }
    }
}

bool StressAutocorrelation::require_data_on_CPU() {
    return false;
}

void StressAutocorrelation::update_data(llint curr_step) {
    if(!_config_info->interaction->has_custom_stress_tensor()) {
        _config_info->interaction->compute_standard_stress_tensor();
    }

    StressTensor stress_tensor = _config_info->interaction->stress_tensor();
    for(uint i = 0; i < stress_tensor.size(); i++) {
        _st_avg[i] += stress_tensor[i];
    }

    // Map used elsewhere in this file:
    // [0]=xx, [1]=yy, [2]=zz, [3]=xy, [4]=xz, [5]=yz
    const double xx = stress_tensor[0];
    const double yy = stress_tensor[1];
    const double zz = stress_tensor[2];
    const double xy = stress_tensor[3];
    const double xz = stress_tensor[4];
    const double yz = stress_tensor[5];

    if(_dump_tensor_cyl && _tensor_cyl_out.is_open()) {
        // Subset = all particles
        const auto &particles = CONFIG_INFO->particles();

        StressTensor st_cyl = _config_info->interaction->stress_tensor_subset_cylindrical(
            particles,
            _cyl_origin,
            _config_info->box->V(),
            true,   // use_min_image
            false   // include_kinetic (match your other dumps)
        );

        const double rr = st_cyl[0];
        const double tt = st_cyl[1]; // hoop
        const double zz2 = st_cyl[2];
        const double rt = st_cyl[3];
        const double rz = st_cyl[4];
        const double tz = st_cyl[5];

        _tensor_cyl_out << curr_step << "  "
                        << rr << " " << rt << " " << rz << "  "
                        << rt << " " << tt << " " << tz << "  "
                        << rz << " " << tz << " " << zz2
                        << "\n";
        _tensor_cyl_out.flush();
    }

    // =========================
    // NEW: dump full 3x3 tensor each update (append-only)
    // =========================
    if(_dump_tensor && _tensor_out.is_open()) {
        _tensor_out << curr_step << "  "
                    << xx << " " << xy << " " << xz << "  "
                    << xy << " " << yy << " " << yz << "  "
                    << xz << " " << yz << " " << zz
                    << "\n";
        _tensor_out.flush(); // optional; remove if you want less IO overhead
    }

    const double P = (xx + yy + zz) / 3.0;

    _sigma_P->add_value(P);
    _sigma_xx->add_value(xx);
    _sigma_yy->add_value(yy);
    _sigma_zz->add_value(zz);
    _sigma_xy->add_value(xy);
    _sigma_yz->add_value(yz);
    _sigma_xz->add_value(xz);

    _N_xy->add_value(xx - yy);
    _N_xz->add_value(xx - zz);
    _N_yz->add_value(yy - zz);

    _times_updated++;
}

// get_output_string unchanged (still prints the ACF-based moduli)
std::string StressAutocorrelation::get_output_string(llint curr_step) {
    std::stringstream ss;

    ss << "# t = " << curr_step << " st_avg =";
    for(uint i = 0; i < _st_avg.size(); i++) {
        ss  << " " << _st_avg[i] / _times_updated;
    }
    ss << std::endl;

    std::vector<double> times;
    _sigma_xy->get_times(_delta_t, times);

    std::vector<double> acf_sigma_xy, acf_sigma_yz, acf_sigma_zx, acf_N_xy, acf_N_xz, acf_N_yz, acf_sigma_xx, acf_sigma_yy, acf_sigma_zz, acf_sigma_P;
    _sigma_P->get_acf(_delta_t, acf_sigma_P);
    _sigma_xx->get_acf(_delta_t, acf_sigma_xx);
    _sigma_yy->get_acf(_delta_t, acf_sigma_yy);
    _sigma_zz->get_acf(_delta_t, acf_sigma_zz);
    _sigma_xy->get_acf(_delta_t, acf_sigma_xy);
    _sigma_yz->get_acf(_delta_t, acf_sigma_yz);
    _sigma_xz->get_acf(_delta_t, acf_sigma_zx);

    _N_xy->get_acf(_delta_t, acf_N_xy);
    _N_xz->get_acf(_delta_t, acf_N_xz);
    _N_yz->get_acf(_delta_t, acf_N_yz);

    double V = _config_info->box->V();
    double T = _config_info->temperature();
    for(uint i = 0; i < times.size(); i++) {
        double bulk_Gt = V / T * acf_sigma_P[i];

        ss << times[i];
        if(_anisotropic) {
            double Gxy = V / T * acf_sigma_xy[i];
            double Gyz = V / T * acf_sigma_yz[i];
            double Gzx = V / T * acf_sigma_zx[i];

            double longitudinal_Gxx = V / T * acf_sigma_xx[i];
            double longitudinal_Gyy = V / T * acf_sigma_yy[i];
            double longitudinal_Gzz = V / T * acf_sigma_zz[i];

            ss << Utils::sformat(" %.25e", Gxy)
               << Utils::sformat(" %.25e", Gyz)
               << Utils::sformat(" %.25e", Gzx)
               << Utils::sformat(" %.25e", bulk_Gt)
               << Utils::sformat(" %.25e", longitudinal_Gxx)
               << Utils::sformat(" %.25e", longitudinal_Gyy)
               << Utils::sformat(" %.25e", longitudinal_Gzz)
               << std::endl;
        }
        else {
            double Gt = V / (5. * T) * (acf_sigma_xy[i] + acf_sigma_yz[i] + acf_sigma_zx[i]);
            Gt += V / (30. * T) * (acf_N_xy[i] + acf_N_xz[i] + acf_N_yz[i]);

            double longitudinal_Gt = V / (3. * T) * (acf_sigma_xx[i] + acf_sigma_yy[i] + acf_sigma_zz[i]);

            ss << Utils::sformat(" %.25e", Gt)
               << Utils::sformat(" %.25e", bulk_Gt)
               << Utils::sformat(" %.25e", longitudinal_Gt)
               << std::endl;
        }
    }

    return ss.str();
}