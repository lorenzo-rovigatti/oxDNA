import gc
import multiprocessing as mp
import sys
from multiprocessing import Array, Process

import numpy as np
import oxpy

import ffs_utils as fu


def save_flux_progress(ctrl):
    with ctrl.progress_lock:
        data = {
            "mode": "flux",
            "success_count": int(ctrl.success_count.value),
            "success_times_sum": int(sum(ctrl.success_times)),
            "success_times_sqr_sum": int(sum(ctrl.success_times_sqr)),
            "completed": bool(ctrl.success_count.value >= ctrl.options["desired_success_count"]),
        }
        fu.save_progress_atomic(data)


def _run_one_flux_restart(idx, ctrl, rng, model):
    input_filename = ctrl.options["input"]
    interface = ctrl.interface

    with oxpy.Context():
        my_input = oxpy.InputFile()
        my_input.init_from_filename(input_filename)
        my_input["print_conf_interval"] = "1e20"
        my_input["print_energy_every"] = "1e20"
        for key in list(my_input.keys()):
            if "data_output" in key:
                del my_input[key]
        for key in ["log_file", "trajectory_file", "lastconf_file", "energy_file"]:
            if key in my_input.keys():
                my_input[key] = "/dev/null"
        my_input["seed"] = str(rng.integers(0, 2**31))
        fu.log("seed: " + my_input["seed"], process=idx)
        my_input["data_output_1"] = ctrl.options["conf_obs_string"]
        my_input["refresh_vel"] = "true"
        my_input["restart_step_counter"] = "true"

        manager = oxpy.OxpyManager(my_input)
        config_info = manager.config_info()
        conf_obs = config_info.get_observable_by_id("success_conf")
        pairs = model.build_pairs(config_info)

        current_state, _ = model.state_ops(config_info, pairs)
        A_recent = current_state == ctrl.options["state_A"]
        dt = 0

        while ctrl.success_count.value < ctrl.options["desired_success_count"]:
            manager.run(1)
            dt += 1
            new_state, new_ops = model.state_ops(config_info, pairs)
            new_ops_string = ", ".join([f"{k} == {v}" for k, v in new_ops.items()])
            if manager.current_step % 1_000_000 == 0:
                fu.log(new_ops_string, process=idx, step=manager.current_step)
            if new_state != current_state:
                fu.log(f"crossing {current_state} --> {new_state} ({new_ops_string})", process=idx, step=manager.current_step, if_verbose=True)
                if new_state == ctrl.options["state_A"]:
                    A_recent = True
                elif interface.crossed_or_overshot(current_state, new_state):
                    if A_recent:
                        with ctrl.success_count.get_lock(), ctrl.success_times.get_lock(), ctrl.success_times_sqr.get_lock():
                            success_idx = ctrl.success_count.value
                            conf = conf_obs.get_output_string(dt)
                            with open(f"success_{success_idx}.dat", "w") as f:
                                f.write(conf + "\n")
                            ctrl.success_count.value += 1
                            ctrl.success_times[idx] += dt
                            ctrl.success_times_sqr[idx] += dt**2
                            fu.log(f"interface reached going forwards ({ctrl.success_count.value} successes)", process=idx, step=manager.current_step)
                        meta = model.success_metadata(config_info, pairs, success_idx, manager.current_step, new_state)
                        if meta is not None:
                            meta["process"] = idx
                            meta["seed"] = int(my_input["seed"])
                            with ctrl.metadata_lock:
                                fu.append_jsonl(ctrl.options["success_metadata_file"], meta)
                        save_flux_progress(ctrl)
                        return "success"
                    fu.log("interface reached without having gone back to A", process=idx, step=manager.current_step)
                elif new_state > interface.Q_right:
                    fu.log(f"crossed to a later state (Q == {new_state}), stopping and restarting", process=idx, step=manager.current_step)
                    return "restart"
                current_state = new_state
        return "done"


def simulate(idx, ctrl):
    fu.log.init(ctrl.options)
    model = fu.build_model(ctrl.options)
    rng = np.random.default_rng([ctrl.options["initial_seed"], idx])
    while ctrl.success_count.value < ctrl.options["desired_success_count"]:
        outcome = _run_one_flux_restart(idx, ctrl, rng, model)
        gc.collect()
        if outcome == "done":
            return


if __name__ == "__main__":
    mp.set_start_method("spawn", force=True)
    if len(sys.argv) < 2:
        print(f"Usage is {sys.argv[0]} input")
        exit(1)
    ctrl = fu.FFSController(sys.argv[1])
    prog = fu.load_or_init_flux_progress(ctrl.options["desired_success_count"])
    if prog["completed"]:
        sys.exit(0)
    ctrl.success_count.value = int(prog["success_count"])
    ctrl.success_times = Array("q", ctrl.options["ncpus"])
    ctrl.success_times_sqr = Array("q", ctrl.options["ncpus"])
    ctrl.success_times[0] = int(prog["success_times_sum"])
    ctrl.success_times_sqr[0] = int(prog["success_times_sqr_sum"])
    fu.log(f"FLUX resume status: success_count={ctrl.success_count.value}, completed={prog['completed']}")
    processes = [Process(target=simulate, args=(i, ctrl)) for i in range(ctrl.options["ncpus"])]
    tp = Process(target=fu.timer, args=(ctrl.success_count, ctrl.options))
    tp.start()
    fu.log(f"Main: Starting {ctrl.options['ncpus']} processes...")
    for p in processes:
        p.start()
    fu.log("Main: waiting for processes to finish")
    failed = False
    for p in processes:
        p.join()
        if p.exitcode != 0:
            failed = True
            fu.log(f"ERROR: worker {p.name} exited with code {p.exitcode}")
    fu.log("Main: Terminating timer")
    tp.kill()
    save_flux_progress(ctrl)
    if failed:
        fu.log("ERROR: one or more FLUX workers failed.")
        sys.exit(1)
    if ctrl.success_count.value < ctrl.options["desired_success_count"]:
        fu.log(f"ERROR: FLUX ended with {ctrl.success_count.value} successes, but desired_success_count={ctrl.options['desired_success_count']}.")
        sys.exit(1)
    avg_timesteps = sum(ctrl.success_times) / ctrl.success_count.value
    avg_timesteps_sqr = sum(ctrl.success_times_sqr) / ctrl.success_count.value
    std_dev = max(0.0, avg_timesteps_sqr - avg_timesteps**2) ** 0.5
    fu.log("-------")
    fu.log("SUMMARY")
    fu.log("-------")
    fu.log(f"Average number of timesteps taken to reach a success (aka 1 / flux): {avg_timesteps}")
    fu.log(f"Average number of timesteps standard deviation: {std_dev}")
    fu.log(f"Initial flux: {1 / avg_timesteps}")
