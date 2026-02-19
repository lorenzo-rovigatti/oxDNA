#ifndef STRESSINBOX_H_
#define STRESSINBOX_H_

#include "BaseObservable.h"

#include <fstream>
#include <string>
#include <vector>

class StressInBox : public BaseObservable {
public:
    StressInBox();
    virtual ~StressInBox();

    void get_settings(input_file &my_inp, input_file &sim_inp) override;
    void init() override;

    // We need positions/forces available on CPU for this observable
    bool require_data_on_CPU() override;

    void update_data(llint curr_step) override;
    std::string get_output_string(llint curr_step) override;

private:
    // Box definition
    LR_vector _center = LR_vector(0., 0., 0.);
    double _L     = 10.0;   // side length
    double _halfL = 5.0;    // computed
    double _Vbox  = 1000.0; // computed

    // Output tensor components (per update)
    double xx = 0.0, yy = 0.0, zz = 0.0;
    double xy = 0.0, xz = 0.0, yz = 0.0;
    double yx = 0.0, zx = 0.0, zy = 0.0;

    // Diagnostics
    int    n_in = 0;
    double Fsum_abs = 0.0;  // sum |F|
    double Fmax_abs = 0.0;  // max |F|
    int    nFpos = 0;       // count |F| > 0

    // Cancellation-proof diagnostics
    LR_vector Fsum_vec = LR_vector(0.,0.,0.); // sum F vector
    double dr_abs_sum = 0.0;                  // sum |dr|
    double dr_abs_max = 0.0;                  // max |dr|
    double virial_abs_sum = 0.0;              // sum |dr_α F_β| over αβ

    // Debug controls
    bool _debug = false;
    llint _debug_every = 1000;
    int _debug_samples = 5;
    std::string _debug_filename = "stress_in_box_debug.dat";
    std::ofstream _dbg;

private:
    inline bool _inside_box(const LR_vector &dr) const {
        return (dr.x >= -_halfL && dr.x <= _halfL &&
                dr.y >= -_halfL && dr.y <= _halfL &&
                dr.z >= -_halfL && dr.z <= _halfL);
    }
};

#endif /* STRESSINBOX_H_ */