/**
 * @file    SimManager.h
 * @date    03/set/2010
 * @author  lorenzo
 *
 *
 */

#ifndef SIMMANAGER_H_
#define SIMMANAGER_H_

#include <cstring>
#include <ctime>

#include "../defs.h"
#include "../Backends/SimBackend.h"
#include "../Utilities/time_scales/time_scales.h"

struct double4;
struct float4;

/**
 * @brief Manages a simulation, be it MC, MD, on GPU or on CPU.
 *
 * @verbatim
[print_input = <bool> (make oxDNA write the input key=value pairs used by the simulation in a file named input.pid, with pid being the oxDNA pid. Defaults to False.)]
conf_file = <string> (path to the starting configuration)

steps = <int> (length of the simulation, in time steps)
[equilibration_steps = <int> (number of equilibration steps. During equilibration, oxDNA does not generate any output. Defaults to 0)]
[restart_step_counter = <bool> (if True oxDNA will reset the step counter to 0, otherwise it will start from the step counter found in the initial configuration. Defaults to False.)]

time_scale = linear/log_lin (a linear time_scale will make oxDNA print linearly-spaced configurations. a log_lin will make it print linearly-spaced cycles of logarithmically-spaced configurations.)
print_conf_interval = <int> (if the time scale is linear, this is the number of time steps between the outputing of configurations, otherwise this is just the first point of the logarithmic part of the log_lin time scale)
print_conf_ppc = <int> (mandatory only if time_scale == log_line. This is the number of printed configurations in a single logarithmic cycle.)
[print_energy_every = <int> (number of time steps between the outputing of the energy (and of the other default observables such as acceptance ratios in Monte Carlo simulations). Defaults to 0.)]
@endverbatim
 */
class SimManager {
protected:
	std::shared_ptr<SimBackend> _backend;
	input_file _input;
	time_scale _time_scale_manager;
	int _time_scale;
	llint _steps, _equilibration_steps;
	llint _steps_run;
	int _seed;
	int _print_energy_every;
	int _pid;
	int _print_input;
	int _fix_diffusion_every;

public:
	SimManager(input_file input);
	virtual ~SimManager();

	static bool stop;
	static bool started;
	virtual void load_options();
	virtual void init();
	virtual void run();
};

void gbl_termhandler(int);

#endif /* SIMMANAGER_H_ */
