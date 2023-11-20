#!/usr/bin/env python3

import sys
import threading
import queue
import os
import subprocess as sp
import math
import difflib
import distutils

from multiprocessing import Lock

SUFFIX_INPUT = "_input"
SUFFIX_COMPARE = "_compare"
SUFFIX_LOG = "_log.dat"

def get_log_name(level):
    return "%s%s" % (level, SUFFIX_LOG)

# static class
class Logger():
    debug_level = 1
    DEBUG = 0
    INFO = 1
    WARNING = 2
    CRITICAL = 3
    RESULTS = 4
    log_lock = Lock()

    messages = ("DEBUG", "INFO", "WARNING", "CRITICAL", "RESULTS")

    @staticmethod
    def log(msg, level, prepend=""):
        if level < Logger.debug_level: 
            return

        Logger.log_lock.acquire()
        print("%s%s: %s" % (prepend, Logger.messages[level], msg), file=sys.stderr)
        Logger.log_lock.release()


class Runner(threading.Thread):
    queue = queue.Queue(1)
    
    def __init__(self, tid):
        threading.Thread.__init__(self)
        self.tid = tid
        self.base_folder = os.getcwd()
        
    def run(self):
        while True:
            details = Runner.queue.get(True)
            
            system = details["system"]
            folder = system.folder
            log_file = get_log_name(details["level"])
            if details["level"] == "oxpy":
                to_execute = f"{details['executable']} {system.input_name}.py"
            else:
                to_execute = f"{details['executable']} {system.input_name} log_file={log_file} no_stdout_energy=0" 
            
            try:
                p = sp.Popen(to_execute, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, cwd=folder, universal_newlines=True)
                p.wait()
                system.simulation_done(p)
                
            except Exception as e:
                Logger.log(e, Logger.WARNING)
            
            Runner.queue.task_done()
            
            
class BaseTest(object):
    def __init__(self, folder, log_prefix, parameters):
        self.log_prefix = log_prefix
        self.folder = folder
        self.parameters = parameters
        self.error = False
        self.parse_parameters()
        
    def parse_parameters(self):
        Logger.log("%s %s's validate_parameters() not implemented" % (self.log_prefix, type(self)), Logger.CRITICAL)
        sys.exit(1)
    
    def test(self):
        Logger.log("%s %s's test() not implemented" % (self.log_prefix, type(self)), Logger.CRITICAL)
        sys.exit(1)
    
    def generate_compare_line(self):
        Logger.log("%s %s's generate_compare_line() not implemented" % (self.log_prefix, type(self)), Logger.CRITICAL)
        sys.exit(1)
    
    
class FileExists(BaseTest):
    def __init__(self, folder, log_prefix, parameters):
        BaseTest.__init__(self, folder, log_prefix, parameters)
        
    def parse_parameters(self):
        if len(self.parameters) < 2:
            Logger.log("%s FileExists expects a single parameter: the name of the file whose existence should be checked" % self.log_prefix, Logger.WARNING)
            self.error = True
        else:
            self.target_file = os.path.join(self.folder, self.parameters[1])
            if len(self.parameters) > 2:
                self.check_if_empty = distutils.util.strtobool(self.parameters[2])
            else:
                self.check_if_empty = True
    
    def test(self):
        if os.path.exists(self.target_file):
            if self.check_if_empty:
                return os.stat(self.target_file).st_size != 0
            return True
        
        return False


class DiffFiles(BaseTest):
    def __init__(self, folder, log_prefix, parameters):
        BaseTest.__init__(self, folder, log_prefix, parameters)
        
    def parse_parameters(self):
        if len(self.parameters) != 3:
            Logger.log("%s DiffFiles expects 2 parameters: the names of the reference and data files, in this order" % self.log_prefix, Logger.WARNING)
            self.error = True
        else:
            self.reference_file = os.path.join(self.folder, self.parameters[1])
            self.data_file = os.path.join(self.folder, self.parameters[2])
    
    def test(self):
        if not self.error:
            with open(self.reference_file) as ref_file:
                ref_data = ref_file.readlines()
            with open(self.data_file) as data_file:
                data = data_file.readlines()
    
        return len(list(difflib.unified_diff(ref_data, data))) == 0
    
    
class ColumnAverage(BaseTest):
    def __init__(self, folder, log_prefix, parameters):
        BaseTest.__init__(self, folder, log_prefix, parameters)
            
    def parse_parameters(self):
        if len(self.parameters) != 5:
            Logger.log("%s ColumnAverage expects 4 parameters: the reference file, the column index, the reference value and its associated error" % self.log_prefix, Logger.WARNING)
            self.error = True
            
        if not self.error:
            self.filename = os.path.abspath(os.path.join(self.folder, self.parameters[1]))
            self.column = int(self.parameters[2]) - 1
            self.true_avg = float(self.parameters[3])
            self.tolerance = float(self.parameters[4])
            
    def _get_average(self):
        if not os.path.exists(self.filename):
            Logger.log("%s ColumnAverage expects the reference file '%s' to exist and be readable" % (self.log_prefix, self.filename), Logger.WARNING)
            self.error = True
            return 0., 0.
        
        avg = 0.
        avg2 = 0.
        c = 0
        
        f = open(self.filename, "r")
        for l in f.readlines():
            v = float(l.split()[self.column])
            avg += v
            avg2 += v*v
            c += 1
        f.close()
        
        avg /= c
        avg2 /= c
        sigma = avg2 - avg*avg
        
        return avg, math.sqrt(sigma)
        
    def test(self):
        avg, error = self._get_average()
        
        if self.error: 
            return False
        
        lower_b = self.true_avg - self.tolerance
        higher_b = self.true_avg + self.tolerance
        
        if lower_b <= avg <= higher_b: 
            Logger.log("%s ColumnAverage test on file '%s', column %d passed. Compute value: %f, Reference value: %f (+- %f)" % (self.log_prefix, self.filename, self.column, avg, self.true_avg, self.tolerance), Logger.DEBUG)
            return True
        
        Logger.log("%s ColumnAverage test on file '%s', column %d failed. Compute value: %f, Reference value: %f (+- %f)" % (self.log_prefix, self.filename, self.column, avg, self.true_avg, self.tolerance), Logger.WARNING)
        
        return False
    
    def generate_compare_line(self):
        avg, error = self._get_average()
        
        spl_line = self.parameters[:]
        spl_line[3] = avg
        spl_line[4] = error
        
        return "::".join(str(x) for x in spl_line)
    
    
class Analyser(object):
    def __init__(self, folder, level):
        self.log_prefix = "Analyser '%s':" % folder
        Logger.log("%s initialising" % self.log_prefix, Logger.DEBUG)

        self.folder = folder        
        self.compare_file = os.path.abspath(os.path.join(folder, level + SUFFIX_COMPARE))
        if not os.path.exists(self.compare_file) or not os.path.isfile(self.compare_file):
            Logger.log("%s compare file '%s' not found" % (self.log_prefix, self.compare_file), Logger.CRITICAL)
            sys.exit(1)

        self.parse_compare_file()
        
    def parse_compare_file(self):
        self.tests = []
        cmp = open(self.compare_file, "r")
        for l in cmp.readlines():
            spl = l.strip().split("::")
            if len(spl) > 0:
                # get an instance of the required test class
                try:
                    test_class = globals()[spl[0]]
                except KeyError as e:
                    Logger.log("%s test class '%s' not found" % (self.log_prefix, spl[0]), Logger.WARNING)
                    continue
                
                self.tests.append(test_class(self.folder, self.log_prefix, spl))
                
        cmp.close()
        
    def print_compare_file(self):
        lines = []
        for test in self.tests:
            lines.append(test.generate_compare_line())
            
        cmp = open(self.compare_file, "w")
        cmp.write("\n".join(lines))
        cmp.close()
            
    
    def test(self):
        n_tests = len(self.tests)
        n_passed = 0
        
        for test in self.tests:
            if test.test(): 
                n_passed += 1
        
        return (n_tests, n_passed)


class System(object):
    def __init__(self, folder, level, exec_name):
        self.log_prefix = "System '%s':" % folder
        
        Logger.log("%s initialising" % self.log_prefix, Logger.DEBUG)
        
        self.executable_name = exec_name
        self.folder = os.path.abspath(folder)
        self.level = level
        self.input_name = level + SUFFIX_INPUT
        self.input_file = os.path.abspath(os.path.join(folder, level + SUFFIX_INPUT))
        
        self.analyser = Analyser(folder, level)
        
        self.error = False
        self.n_tests = 0
        self.n_passed = 0
    
    def simulation_done(self, p, do_tests=True):
        self.error = False
        if p.returncode != 0:
            # segfault
            if p.returncode == 139: 
                Logger.log("%s segmentation fault (return code %d)" % (self.log_prefix, p.returncode), Logger.WARNING)
            else: 
                Logger.log("%s %s returned %d" % (self.executable_name, self.log_prefix, p.returncode), Logger.WARNING)
            self.error = True
            
        # check the logfile for errors
        log_file = os.path.join(self.folder, get_log_name(self.level))
        if os.path.exists(log_file):
            f = open(log_file)
            for l in f.readlines():
                if l.startswith("ERROR"): 
                    Logger.log("%s %s error: %s" % (self.executable_name, self.log_prefix, l.strip()), Logger.WARNING)
                    self.error = True
            f.close()
        
        # check the standard output for nans and infs
        for l in p.stdout.readlines():
            # we need byte-like objects and not 'str' 
            tokens = ["nan", "inf"]
            for t in tokens: 
                if t in l: 
                    Logger.log("%s %s generated a '%s': %s" % (self.executable_name, t, l), Logger.WARNING)
                    self.error = True
        
        # we don't run tests if the simulation was not successful. We put this here so that all
        # above messages can be printed independently of each other
        if self.error: 
            return
        
        Logger.log("%s %s run completed and successful" % (self.log_prefix, self.executable_name), Logger.DEBUG)
        
        if do_tests:
            (n_tests, n_passed) = self.analyser.test()
            Logger.log("%s\n\tPassed/Total: %d/%d" % (self.log_prefix, n_passed, n_tests), Logger.RESULTS)
            
            self.n_tests = n_tests
            self.n_passed = n_passed
        
    
class TestManager(object):
    def __init__(self, list_file, executable, level, threads=1):
        self.log_prefix = "TestManager:"
        
        self.executable = executable
        self.executable_name = os.path.basename(self.executable)
        self.level = level
        self.systems = []
        self.threads = threads
            
        input_name = level + SUFFIX_INPUT
        if level == "oxpy":
            input_name += ".py"
        
        f = open(list_file, "r")
        for l in f.readlines():
            d = l.strip()
            if len(d) == 0 or d[0] == '#':
                continue
            input_path = os.path.join(d, input_name)
            if os.path.exists(input_path) and os.path.isfile(input_path):
                new_system = System(d, level, self.executable_name) 
                self.systems.append(new_system)
                
        f.close()
                
    def launch(self):
        for i in range(self.threads):
            runner = Runner(i)
            runner.daemon = True
            runner.start()
        
        for system in self.systems:
            Logger.log("%s starting %s in '%s'" % (self.executable_name, self.log_prefix, system.folder), Logger.DEBUG)
            
            details = {
                       "system" : system,
                       "executable" : self.executable,
                       "level" : self.level
                       }
            
            Runner.queue.put(details, block=True)
            
        Runner.queue.join()
        
    def finalise(self):
        n_passed = 0
        n_tests = 0
        n_errors = 0
        for system in self.systems:
            n_passed += system.n_passed
            n_tests += system.n_tests
            if system.error:
                n_errors += 1
                
        if n_errors == 1:
            Logger.log("The analysis in %d folder failed" % n_errors, Logger.CRITICAL, "\n")
        elif n_errors > 1:
            Logger.log("The analysis in %d folders failed" % n_errors, Logger.CRITICAL, "\n")
            
        Logger.log("Summary for level '%s'\n\tPassed/Total: %d/%d\n" % (self.level, n_passed, n_tests), Logger.RESULTS, "\n")
        
        if n_tests != n_passed or n_errors > 0:
            Logger.log("Not all tests have passed successfully", Logger.CRITICAL)
    
    
if __name__ == '__main__':
    import argparse
    
    def print_version():
        print("oxDNA Test Suite v 0.1")
        print("Copyright (C) 2015 oxDNA")
        print("This is free software; see the source for copying conditions.  There is NO")
        print("warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n")
        exit(1)
        
    parser = argparse.ArgumentParser(description="A simple script that runs tests in selected folders and compare the results with reference data.")
    parser.add_argument("--debug", action="store_true", help="Print debug messages")
    # parser.add_argument("--help", action="store_false", help="Print this message")
    parser.add_argument("--version", action="store_true", help="Print the script version and exit")
    parser.add_argument("--threads", default=1, help="Set the number of concurrent threads that will be launched")
    parser.add_argument("folder_list_file")
    parser.add_argument("executable")
    parser.add_argument("test_level")
    
    args = parser.parse_args()
    
    if args.debug:
        Logger.debug_level = 0
    if args.version:
        print_version()
        exit(0)
        
    if not os.path.exists(args.folder_list_file) or not os.path.isfile(args.folder_list_file):
        Logger.log("List file '%s' does not exist or it is unreadable" % args.folder_list_file, Logger.CRITICAL)
        sys.exit(1)
        
    Logger.log("Running tests for level '%s'" % sys.argv[3], Logger.INFO)
    tm = TestManager(args.folder_list_file, args.executable, args.test_level, threads=args.threads)
    tm.launch()
    tm.finalise()
    