import glob
import json
import math
import os
import sys
import time
from multiprocessing import Lock, Value

import oxpy

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib


class Logger:
    def __init__(self):
        self.log_lock = Lock()
        self.verbose = False
        self.also_stderr = False
        self.initialised = False
        self.log_file = None

    def init(self, options):
        if "log" not in options or "filename" not in options["log"]:
            print("Required key 'log.filename' not found, exiting", file=sys.stderr)
            exit(1)
        self.log_file = open(options["log"]["filename"], "a", buffering=1)
        self.verbose = bool(options["log"].get("verbose", False))
        self.also_stderr = bool(options["log"].get("also_stderr", False))
        self.initialised = True

    def __call__(self, text, process=None, step=None, if_verbose=False):
        if if_verbose and not self.verbose:
            return
        if not self.initialised:
            print("Logger uninitialised, exiting", file=sys.stderr)
            exit(1)
        if step is not None:
            text = f"step {step}, {text}"
        if process is not None:
            text = f"Process {process}, {text}"
        with self.log_lock:
            self.log_file.write(text + "\n")
            self.log_file.flush()
            if self.also_stderr:
                print(text, file=sys.stderr)


log = Logger()


def _require(cond, msg):
    if not cond:
        log(msg)
        exit(1)


def append_jsonl(filename, record):
    with open(filename, "a") as f:
        f.write(json.dumps(record) + "\n")
        f.flush()
        os.fsync(f.fileno())


PROGRESS_FILE = "progress.json"


def load_progress():
    if not os.path.exists(PROGRESS_FILE):
        return {}
    with open(PROGRESS_FILE, "r") as f:
        return json.load(f)


def save_progress_atomic(data):
    tmp = PROGRESS_FILE + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, PROGRESS_FILE)


def count_existing_success_files():
    max_idx = -1
    for fn in glob.glob("success_*.dat"):
        try:
            stem = os.path.splitext(os.path.basename(fn))[0]
            idx = int(stem.split("_")[1])
            max_idx = max(max_idx, idx)
        except Exception:
            pass
    return max_idx + 1


def load_or_init_flux_progress(desired_success_count):
    prog = load_progress()
    existing_successes = count_existing_success_files()
    success_count = max(int(prog.get("success_count", 0)), existing_successes)
    completed = bool(prog.get("completed", False)) or success_count >= desired_success_count
    return {
        "mode": "flux",
        "success_count": success_count,
        "success_times_sum": int(prog.get("success_times_sum", 0)),
        "success_times_sqr_sum": int(prog.get("success_times_sqr_sum", 0)),
        "completed": completed,
    }


def load_or_init_shoot_progress(total_simulations, n_starting_confs):
    prog = load_progress()
    existing_successes = count_existing_success_files()
    attempt_count = int(prog.get("attempt_count", 0))
    success_count = max(int(prog.get("success_count", 0)), existing_successes)
    completed = bool(prog.get("completed", False)) or attempt_count >= total_simulations
    success_from = prog.get("success_from", [0] * n_starting_confs)
    attempt_from = prog.get("attempt_from", [0] * n_starting_confs)
    if len(success_from) != n_starting_confs:
        success_from = [0] * n_starting_confs
    if len(attempt_from) != n_starting_confs:
        attempt_from = [0] * n_starting_confs
    return {
        "mode": "shoot",
        "attempt_count": attempt_count,
        "success_count": success_count,
        "success_from": success_from,
        "attempt_from": attempt_from,
        "completed": completed,
    }


class Interface:
    def __init__(self, Q_left, Q_right):
        if Q_right - Q_left != 1:
            log(
                f"Interfaces separating states further away than 1 are not supported "
                f"(Q_left == {Q_left}, Q_right == {Q_right})"
            )
            exit(1)
        self.Q_left = Q_left
        self.Q_right = Q_right

    def crossed(self, Q1, Q2):
        return Q1 == self.Q_left and Q2 == self.Q_right

    def crossed_or_overshot(self, Q1, Q2):
        return Q1 == self.Q_left and Q2 >= self.Q_right


class Options(dict):
    default_options = {
        "distance_final_state": 0,
        "bond_threshold": -0.10,
        "use_bonds": True,
        "success_metadata_file": "success_metadata.jsonl",
        "conf_obs_string": """{
  name = /dev/null
  print_every = 1e20
  only_last = true
  col_1 = {
    type = Configuration
    id = success_conf
  }
}""",
    }

    def __init__(self, input_file):
        super().__init__(Options.default_options)
        try:
            with open(input_file, "rb") as f:
                self.update(tomllib.load(f))
        except FileNotFoundError:
            print("Input file not found or not accessible, exiting", file=sys.stderr)
            exit(1)
        except tomllib.TOMLDecodeError:
            print("Input file contains non-valid TOML, exiting", file=sys.stderr)
            exit(1)

    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except KeyError:
            log(f"Option '{key}' not found, exiting")
            exit(1)


class FFSController:
    def __init__(self, input_file):
        self.options = Options(input_file)
        log.init(self.options)
        self.interface = Interface(int(self.options["interface"][0]), int(self.options["interface"][1]))
        self.success_count = Value("i", 0)
        self.progress_lock = Lock()
        self.metadata_lock = Lock()


def timer(success_count, options=None):
    if options is not None and not log.initialised:
        log.init(options)
    timestep = time.asctime(time.localtime())
    log("Timer started at %s" % timestep)
    itime = time.time()
    while True:
        time.sleep(10)
        now = time.time()
        with success_count.get_lock():
            ns = success_count.value
            if ns > 0:
                log("Timer: at %s: successes: %d, time per success: %g (%g sec elapsed)" % (timestep, ns, (now - itime) / ns, now - itime))
            else:
                log("Timer: at %s: no successes yet" % timestep)


class HybridPairOP:
    def __init__(self, options):
        self.options = options
        _require("pairs" in options, "Missing required key 'pairs' in ffs.toml.")
        self.pair_indices = [(int(i), int(j)) for i, j in options["pairs"]]
        _require(len(self.pair_indices) > 0, "At least one pair must be provided.")
        _require("distance_thresholds" in options, "Missing required key 'distance_thresholds' in ffs.toml.")
        self.distance_thresholds = [float(x) for x in options["distance_thresholds"]]
        _require(len(self.distance_thresholds) > 0, "distance_thresholds cannot be empty.")
        _require(all(self.distance_thresholds[i] > self.distance_thresholds[i + 1] for i in range(len(self.distance_thresholds) - 1)), "distance_thresholds must be strictly decreasing.")
        self.distance_final_state = int(options.get("distance_final_state", 0))
        self.use_bonds = bool(options.get("use_bonds", True))
        self.bond_threshold = float(options.get("bond_threshold", -0.10))

    def build_pairs(self, config_info):
        particles = config_info.particles()
        n_particles = len(particles)
        pairs = []
        for i, j in self.pair_indices:
            _require(0 <= i < n_particles and 0 <= j < n_particles, f"Pair ({i}, {j}) out of range for {n_particles} particles.")
            pairs.append((particles[i], particles[j]))
        log(f"HybridPairOP pairs={self.pair_indices}")
        log("Distance thresholds=%s, distance_final_state=%d" % (self.distance_thresholds, self.distance_final_state))
        log(f"use_bonds={self.use_bonds}, bond_threshold={self.bond_threshold}")
        return pairs

    def _pair_distances(self, config_info, pairs):
        box = config_info.box
        return [math.sqrt(box.sqr_min_image_distance(p0, p1)) for p0, p1 in pairs]

    def _hb_energies(self, config_info, pairs):
        interaction = config_info.interaction
        interaction.begin_energy_computation()
        return [interaction.pair_interaction_term(4, p0, p1) for p0, p1 in pairs]

    def _distance_state(self, mindist):
        n = len(self.distance_thresholds)
        state = self.distance_final_state
        for k, threshold in enumerate(self.distance_thresholds):
            if mindist > threshold:
                state = self.distance_final_state - n + k
                break
        return state

    def state_ops(self, config_info, pairs):
        distances = self._pair_distances(config_info, pairs)
        mindist = min(distances)
        energies = []
        bonded_indices = []
        nbonds = 0
        if self.use_bonds:
            energies = self._hb_energies(config_info, pairs)
            bonded_indices = [idx for idx, energy in enumerate(energies) if energy < self.bond_threshold]
            nbonds = len(bonded_indices)
        if self.use_bonds and nbonds > 0:
            state = nbonds
        else:
            state = self._distance_state(mindist)
        bonded_pairs = [self.pair_indices[idx] for idx in bonded_indices]
        return state, {
            "Q": state,
            "mindist": mindist,
            "pair_distances": distances,
            "use_bonds": self.use_bonds,
            "nbonds": nbonds,
            "bonded_pair_indices": bonded_indices,
            "bonded_pairs": bonded_pairs,
            "hb_energies": energies,
        }

    def success_metadata(self, config_info, pairs, success_idx, step, q_reached):
        _, ops = self.state_ops(config_info, pairs)
        return {"success_idx": int(success_idx), "step": int(step), "Q": int(q_reached), "pairs": self.pair_indices, "ops": ops}


def build_model(options):
    return HybridPairOP(options)
