import gc
import glob
import multiprocessing as mp
import sys
from multiprocessing import Array, Process, Queue, Value

import numpy as np
import oxpy

import ffs_utils as fu


def save_shoot_progress(ctrl):
    with ctrl.progress_lock:
        data = {
            "mode": "shoot",
            "attempt_count": int(ctrl.attempt_count.value),
            "success_count": int(ctrl.success_count.value),
            "success_from": list(ctrl.success_from[:]),
            "attempt_from": list(ctrl.attempt_from[:]),
            "completed": bool(ctrl.attempt_count.value >= ctrl.options["total_simulations"]),
        }
        fu.save_progress_atomic(data)


def _run_one_shoot_attempt(attempt_id, idx, ctrl, rng, model):
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
        fu.log(f"seed: {my_input['seed']} (attempt {attempt_id})", process=idx)
        my_input["data_output_1"] = ctrl.options["conf_obs_string"]
        my_input["refresh_vel"] = "true"
        my_input["restart_step_counter"] = "true"
        conf_index = rng.integers(0, len(ctrl.starting_confs))
        my_input["conf_file"] = ctrl.starting_confs[conf_index]
        with ctrl.attempt_count.get_lock():
            ctrl.attempt_count.value += 1
            ctrl.attempt_from[conf_index] += 1
        save_shoot_progress(ctrl)
        manager = oxpy.OxpyManager(my_input)
        config_info = manager.config_info()
        conf_obs = config_info.get_observable_by_id("success_conf")
        pairs = model.build_pairs(config_info)
        current_state, _ = model.state_ops(config_info, pairs)
        while True:
            manager.run(1)
            new_state, new_ops = model.state_ops(config_info, pairs)
            new_ops_string = ", ".join([f"{k} == {v}" for k, v in new_ops.items()])
            if manager.current_step % 1_000_000 == 0:
                fu.log(new_ops_string, process=idx, step=manager.current_step)
            if new_state != current_state:
                fu.log(f"crossing {current_state} --> {new_state} ({new_ops_string})", process=idx, step=manager.current_step, if_verbose=True)
                if new_state == ctrl.options["state_A"]:
                    fu.log("failure (back to Q == Q_A)", process=idx, step=manager.current_step)
                    return "failure"
                elif interface.crossed_or_overshot(current_state, new_state):
                    with ctrl.success_count.get_lock():
                        success_idx = ctrl.success_count.value
                        conf = conf_obs.get_output_string(manager.current_step)
                        with open(f"success_{success_idx}.dat", "w") as f:
                            f.write(conf + "\n")
                        ctrl.success_count.value += 1
                        ctrl.success_from[conf_index] += 1
                    meta = model.success_metadata(config_info, pairs, success_idx, manager.current_step, new_state)
                    if meta is not None:
                        meta["process"] = idx
                        meta["seed"] = int(my_input["seed"])
                        meta["starting_conf"] = ctrl.starting_confs[conf_index]
                        with ctrl.metadata_lock:
                            fu.append_jsonl(ctrl.options["success_metadata_file"], meta)
                    save_shoot_progress(ctrl)
                    fu.log(f"interface reached ({ctrl.success_count.value} successes)", process=idx, step=manager.current_step)
                    return "success"
                current_state = new_state


def runner(idx, q, ctrl):
    fu.log.init(ctrl.options)
    rng = np.random.default_rng([ctrl.options["initial_seed"], idx])
    model = fu.build_model(ctrl.options)
    while True:
        if q.empty():
            return
        attempt_id = q.get()
        _run_one_shoot_attempt(attempt_id, idx, ctrl, rng, model)
        gc.collect()


if __name__ == "__main__":
    mp.set_start_method("spawn", force=True)
    if len(sys.argv) < 2:
        print(f"Usage is {sys.argv[0]} input")
        exit(1)
    ctrl = fu.FFSController(sys.argv[1])
    ctrl.starting_confs = glob.glob(ctrl.options["starting_conf_pattern"])
    if len(ctrl.starting_confs) < 1:
        print("0 starting configurations! aborting", file=sys.stderr)
        sys.exit(1)
    prog = fu.load_or_init_shoot_progress(ctrl.options["total_simulations"], len(ctrl.starting_confs))
    if prog["completed"]:
        sys.exit(0)
    fu.log("Main: Found %d configurations with pattern: %s" % (len(ctrl.starting_confs), ctrl.options["starting_conf_pattern"]))
    ctrl.attempt_count = Value("i", int(prog["attempt_count"]))
    ctrl.success_count.value = int(prog["success_count"])
    ctrl.success_from = Array("i", len(ctrl.starting_confs))
    ctrl.attempt_from = Array("i", len(ctrl.starting_confs))
    for i, v in enumerate(prog["success_from"]):
        ctrl.success_from[i] = int(v)
    for i, v in enumerate(prog["attempt_from"]):
        ctrl.attempt_from[i] = int(v)
    fu.log(f"SHOOT resume status: attempt_count={ctrl.attempt_count.value}, success_count={ctrl.success_count.value}, completed={prog['completed']}")
    q = Queue()
    for i in range(ctrl.attempt_count.value, ctrl.options["total_simulations"]):
        q.put(i)
    processes = [Process(target=runner, args=(i, q, ctrl)) for i in range(ctrl.options["ncpus"])]
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
    save_shoot_progress(ctrl)
    if failed:
        fu.log("ERROR: one or more SHOOT workers failed.")
        sys.exit(1)
    nsuccesses = ctrl.success_count.value
    fu.log("-------")
    fu.log("SUMMARY")
    fu.log("-------")
    fu.log("nsuccesses: %d nattempts: %d success_prob: %g" % (nsuccesses, ctrl.attempt_count.value, nsuccesses / ctrl.attempt_count.value))
    with open("conf_statistics.log", "w") as out:
        out.write("# log of successes probabilities from each starting conf\n")
        out.write("# idx nsuccesses nattempts prob\n")
        for k, v in enumerate(ctrl.success_from):
            txt = "%3d    %3d    %3d   " % (k, v, ctrl.attempt_from[k])
            if ctrl.attempt_from[k] > 0:
                txt += "%g" % (v / ctrl.attempt_from[k])
            else:
                txt += "NA"
            out.write(txt + "\n")
