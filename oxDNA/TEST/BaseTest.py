import inspect
import logging
import os
import stat
import string
import subprocess
import sys

# Logging
logger = logging.getLogger('BaseTest')
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)

# Safely load oxDNA UTILS
try:
    sys.path.append('../UTILS')
    #Necessary for when Testers are imported directly from the UserTest subfolder
    sys.path.append('../../UTILS')
    import readers
except (ImportError):
    logger.error("oxDNA UTILS directory not found. Use at your own risk. Some functionality will be broken.")


class TestSuite:
    """
    Runs a list of instantiated Tester objects.
    Output result in human readable format.
    """
    def __init__(self, runner_instance_list=[], tester_class_list=[], verbose=False):
        """
        Initialize the TestSuite
        @ Input checking on runner_instance)list and tester_class_list
        @ Create testers pair-wise from runners and testers
        Note that initialized testers can be added manually with TestSuit.add_test(Tester)
        """
        # Logging
        self.verbose = bool(verbose)
        if self.verbose:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.WARNING)

        # Instantiatzed test list
        # There will be a test for each tester_class/runner_instance combination
        self.tester_instance_list = []
        self.tester_instance_list = self.append_tester_instance_list(runner_instance_list, tester_class_list)

    def __iter__(self):
        return self

    def __len__(self):
        return len(self.tester_instance_list)

    def __str__(self):
        return "TestSuite with %s UserTest" % self.__len__()

    def add_test(self, obj_tester):
        """
        Add a test to the tester_instance_list
        """
        if isinstance(obj_tester,Tester):
            self.tester_instance_list.append(obj_tester)

    def append_tester_instance_list(self, runner_instance_list, tester_class_list):
        """
        Return tester_instance_list
        This is a list of all instantiated tests to run.
        """
        for obj_runner in runner_instance_list:
            if isinstance(obj_runner,Runner):
                self.runner_instance_list.append(obj_runner)

        for class_tester in tester_class_list:
            if inspect.isclass(class_tester):
                self.tester_class_list.append(class_tester)

        for tester_class in tester_class_list:
            for runner_instance in runner_instance_list:
                obj_test = tester_class(runner_instance)
                self.add_test(obj_test)
        return self.tester_instance_list

    def next(self):
        """
        Return the next untested (Tester.test_pass=None) TesterInstance in self.tester_instance_list
        """
        for obj_test in self.tester_instance_list:
            if obj_test.test_pass == None:
                return obj_test
        raise StopIteration

    def result(self):
        """
        Print human-readable TestSuite results to STDOUT
        """
        self.N_failed = 0
        self.N_total = len(self.tester_instance_list)
        last_name = ''
        for obj_test in self.tester_instance_list:
            curr_name = obj_test.obj_runner.run_name
            if curr_name != last_name:
                print "Runner: %s" % (obj_test.obj_runner.run_name)
                last_name = curr_name

            # Only print results if fail or verbose
            if obj_test.test_pass and self.verbose==False:
                continue
            elif not obj_test.test_pass:
                self.N_failed += 1

            print "\tTester: %s" % (obj_test)
            print "\t\t%s, %s" % (obj_test.test_result())
            if self.verbose:
                print "\t\t%s" % (obj_test.test_results_dict)

        print "\tFailed/Total: %i/%i" % (self.N_failed, self.N_total)

    def run(self):
        """
        Run all test in the tester_instance_list
        """
        if not self.tester_instance_list:
            logger.critical("No test added. Add test with self.add_test()")
            return False

        # Benefit: Does not re-run test that have already been run
        try:
            while True: self.next().test()
        except StopIteration:
            pass
        #for obj_test in self.tester_instance_list:
        #    obj_test.test()
        return True

class Tester:
    """
    Abstract Base class for Test classes.
    Keeps track of everything so that the end-user only has to inherent this class and write test in do_test()

    USAGE:
    Example
    class TEST_DEMO(Tester):
        def self.do_test(self):
            <User writes test here>

    INPUT:
        obj_runner --- Runner object. Can be either run or unrun
        output_files_path --- Tester generated output will go here. Output specified in the oxDNA INPUT file will not be affected.
        expected_results_path --- Expected results, for comparison will go here.
            Format: Key = Value
            Comments (#) are ignored

    USER VARIABLES:
    self.test_pass (default False) (REQUIRED!!!)
        User must set it to true under their conditions.
    self.test_pass_info --- Extra human readable test info.
    self.test_results_dict --- Dictionary of test results. Written in a file named after the TesterClass in output_results_path
    self.expected_results_dict --- Dictionary of expected results from file (named after Tester Class) in expected_results_path
        Note: Everything is a string by default. Manual cast for numbers.

    ADDITIONAL USER NOTES:
    - The Runner already has trajectory information loaded in Runner.run_trajectory
    - See Runner documentation for Runner USER VARIABLES
    """
    def __init__(self, obj_runner, output_files_path='output_files/', expected_results_path='expected_results/', verbose=False):
        """
        Initialize Tester object
        @ Load the runner with self.__update_runner(obj_runner)
        """
        # Logging
        self.verbose = bool(verbose)
        if self.verbose:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.WARNING)

        # User-variables for do_test()
        self.test_name = None
        self.test_pass = None
        self.test_pass_info = "Tester not inititialized"
        self.test_results_dict= {}
        self.expected_results_dict = {}

        # Update the Runner
        self.obj_runner = None
        self.obj_runner =  self.__update_runner(obj_runner)

        # filename / filepath checking
        if check_path(output_files_path, sys_exit=False):
            self.output_files_path = str(output_files_path)
        else:
            try:
                os.makedirs(output_files_path)
            except OSError:
                logger.critical("Could not create directory (%s)", output_files_path)

        if check_path(expected_results_path, sys_exit=False):
            self.expected_results_path = str(expected_results_path)
        else:
            try:
                os.makedirs(expected_results_path)
            except OSError:
                logger.critical("Could not create directory (%s)", expected_results_path)

    def __str__(self):
        """
        Return the name of the class
        """
        return self.__class__.__name__

    def do_test(self):
        """
        Overload this function to do testing.
        returns True/False/None
        You MUST update self.test_pass and return it
            True --- Test passed
            False --- Test Failed
            None --- Anything else
        """
        raise NotImplementedError

    def test(self):
        """
        Setup and tear down to test.
        Calls self.do_test() to actually do the test.
        do_test must return True/False (self.test_pass)
        @ Read expected results with self.__read_expected_results()
        """
        # Check if the runner was a success
        # If it was, don't bother running again
        if self.obj_runner.run_success != True:
            self.obj_runner.run()
            self.obj_runner = self.__update_runner(self.obj_runner)

        # Print Runner STDERR to file
        self.__print_log()

        # Read expected results from file
        self.__read_expected_results_file()

        # Test defaults to failed
        self.test_pass = False
        if self.obj_runner.run_success == True:
            self.test_pass_info = "Test Failed."
            try:
                self.test_pass = self.do_test()
                self.__print_test_results()
            except Exception, inst:
                self.test_pass_info = "Test Failed (%s). Exception in do_test()" % self.__str__()
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                logger.exception("Exception in Tester (%s) do_test()", self.__str__())
                #logger.error("%s in %s, Line %s", exc_type.__name__, fname, exc_tb.tb_lineno)
                #logger.error("%s", exc_tb.print_tb())
        else:
            pass
        return self.test_pass

    def test_result(self):
        """
        Returns the state of the test.
        Prints test_results to file
        [test_pass, test_pass_info]
        test_pass
            True
            False
            None
        test_pass_info
            Human readable text about the test state
        """
        return self.test_pass, self.test_pass_info

    def __print_log(self):
        """
        Print log file
        Return True/False for successful execution.
        """
        self.log_filename = "log.txt"
        self.log_filepath = os.path.join(self.output_files_path, self.log_filename)
        try:
            file = open(self.log_filepath, 'w')
            for line in self.obj_runner.STDERR:
                file.write(line)
            file.close()
        except:
            logger.warning("Could not print Runner log (%s)", self.log_filepath)
            return False
        return True

    def __print_test_results(self):
        """
        Output self.test_results_dict to file
        """
        self.test_results_filepath = os.path.join(self.output_files_path, self.__str__())
        file = open(self.test_results_filepath, 'w')
        for key, value in self.test_results_dict.items():
            file.write("%s = %s\n" % (key, value))
        file.close()
        return True

    def __read_expected_results_file(self):
        """
        Read the expected_results_file. Update self.expected_results_dict
        default_path: expected_results/<RunnerName>/<TesterName>
        <RunnerName> from Runner.run_name
        <TesterName> from Tester.__str__()
        return True if read OK. Otherwise False.
        """
        self.expected_results_dict = {}
        self.expected_results_filepath = os.path.join(self.expected_results_path, self.__str__())
        if check_file(self.expected_results_filepath, sys_exit=False):
            fin = open(self.expected_results_filepath)
            for line in fin:
                line = line.lstrip()
                if not line.startswith('#'):
                    (key, value) = line.split('=')
                    key = key.strip()
                    value = value.strip()
                    self.expected_results_dict.update({key: value})
            return True
        else:
            logger.warning("Cannot read expected_results (%s)", self.expected_results_filepath)
        return False

    def __update_runner(self, obj_runner):
        """
        Check and load the runner
        Return the obj_runner
        """
        if not isinstance(obj_runner,Runner):
            logger.error("obj_runner (%s) is not a Runner class")
            return None

        if obj_runner.run_success == True:
            self.test_pass = None
            self.test_pass_info = "Runner loaded. Pending Test"
        elif obj_runner.run_success == False:
            self.test_pass = False
            self.test_pass_info = "Runner did not execute successfully"
        elif obj_runner.run_success == None:
            self.test_pass = None
            self.test_pass_info = "Runner has not executed."
        else:
            logger.critical("Problem with load_runner")
            sys.exit()
        return obj_runner


class Runner:
    """
    Manages the running of oxDNA. Actual execution is held until Runner.run()
    INPUT:
        exec_filepath --- Filepath to oxDNA executable
        input_filename --- oxDNA inputfile
        input_path --- Path to oxDNA inputfile. Note that oxDNA is run from this location
    USER VARIABLES:
        self.run_name --- (default self.input_path)
        self.run_success --- True/False/None
            None  - Run has not been attempted.
            True  - Run has completed.
            False - Run has been attempted, but failed.
        self.run_trajectory --- List of base.py systems (NOTE: UTILS.readers dependency)
        self.run_input_options --- Dictionary of selected options from oxDNA INPUT file
    """
    def __init__(self, exec_filepath="../bin/oxDNA", input_filename="INPUT", input_path=".", run_name="", verbose=False):
        """
        Initialize the Runner object
        @ self.run_success = None
        """
        # Logger
        self.verbose = bool(verbose)
        if self.verbose:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.WARNING)

        # Check runner inputs
        if check_path(input_path, sys_exit=False):
            self.input_path = str(input_path)
        else:
            self.run_success = False
            logger.critical("Invalid input_path (%s)", input_path)
            return None

        self.exec_filepath = str(exec_filepath)
        self.exec_abs_filepath = os.path.abspath(self.exec_filepath)
        if not check_file(self.exec_abs_filepath, sys_exit=False):
            self.run_success = False
            logger.critical("Invalid exec_abs_filepath (%s)", self.exec_abs_filepath)
            return None

        self.input_filename = str(input_filename)
        self.input_filepath = os.path.join(self.input_path, self.input_filename)
        if not check_file(self.input_filepath, sys_exit=False):
            self.run_success = False
            logger.critical("Invalid input_filepath (%s)", self.input_filepath)
            return None

        # User variables
        if run_name:
            self.run_name = str(run_name)
        else:
            self.run_name = self.input_filepath
        self.run_success = None
        self.run_trajectory = []

        # Parse inputfile for additional filenames
        self.run_input_options = {
                'topology'       : None,
                'conf_file'      : None,
                'lastconf_file'  : 'last_conf.dat',
                'trajectory_file': None,
                'energy_file'    : None,
        }
        return None

    def __str__(self):
        """
        Return run_name
        """
        return self.run_name

    def run(self):
        """
        Run oxDNA and load trajectory information.
        self.run_trajectory --- Loaded by self.__import_trajectory()
        Parse input file for additional filenames/filepaths
        """
        # Using Popen for now, but running the Runner could easily be batched/distributed.
        p = subprocess.Popen([self.exec_abs_filepath, self.input_filename], cwd=self.input_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = p.communicate()
        self.STDOUT = output[0]
        self.STDERR = output[1]

        # Quick and dirty checking to see if oxDNA ran successfully
        if (string.find(self.STDERR, 'GAME OVER') != -1):
            logger.info("Run Failed (%s)", self.input_filepath)
            self.run_success = False
        if (string.find(self.STDERR, 'END OF THE SIMULATION') > 0):
            self.run_success = True
            self.__parse_inputfile()
            self.__import_trajectory()
        return self.run_success

    def __get_value_from_label(self, label=""):
        """
        Return the value for a given input file label
        Only returns the value for the first instance of the label

        label --- Label in input file
            Ex: label=topology
                Returns the topology filename
        """
        label_value = None
        fin = open(self.input_filepath)
        for line in fin:
            line = line.lstrip()
            if not line.startswith('#'):
                if label in line:
                    label_value = line.split('=')[1].replace(' ','').replace('\n','')
                    break
        return label_value

    def __import_trajectory(self):
        """
        Import system trajectory into Runner object
            self.run_trajectory --- List of systems.
        """
        conffile = self.run_input_options['trajectory_file']
        topologyfile = self.run_input_options['topology']
        try:
            myreader = readers.LorenzoReader(conffile,topologyfile)
            mysystem = myreader.get_system()
        except (NameError):
            logger.error("Trajectory cannot be read (%s)", self.input_filepath)
            return False

        self.trajectory = []
        while mysystem != False:
            self.run_trajectory.append(mysystem)
            mysystem = myreader.get_system()
        return True

    def __parse_inputfile(self):
        """
        Load selected options from inputfile into self.run_input_options
        No return / Called from __init__()
        """
        for key, value in self.run_input_options.items():
            value = self.__get_value_from_label(key)
            # Skip the value if the key is not found in the inputfile
            if not value:
                logger.info("%s not found in oxDNA INPUT file", key)
                continue
            filepath = os.path.join(self.input_path, value)
            check_file(filepath)
            self.run_input_options[key] = filepath


# Utilities...

def check_path(base_path=".", sys_exit=True):
    """
    Check that the base_path exists.
    If not, critical exit.
    """
    if os.path.isdir(base_path):
        return True
    else:
        if sys_exit == True:
            logger.critical("Invalid path (%s)", base_path)
            sys.exit()
        else:
            return False

def check_file(filepath, sys_exit=True):
    """
    Check that the file exists at base_path/filepath
    If not, critical exit.
    Return True/False
    True if file exists and is not empty
    False otherwise

    sys_exit --- True/False
        Turns on/off printing + exiting
    """
    try:
        if os.stat(filepath)[stat.ST_SIZE]==0:
            if sys_exit == True:
                logger.critical("Empty file (%s)", filepath)
                sys.exit()
            else:
                return False
    except (OSError):
        if sys_exit == True:
            logger.critical("File does not exist (%s)", filepath)
            sys.exit()
        else:
            return False
    return True

def load_file(filename):
    """
    Return the contents of a file as a list.
    Each line in one element of the list.
    Used to load RUNNER_LIST and TEST_LIST
        RUNNER_LIST, each line (element in list) is a directory to walk through
        TEST_LIST, each line (element in list) is a Tester class to use
    """
    check_file(filename)
    fin = open(filename)
    dir_list = []
    for line in fin:
        line = line.lstrip()
        if not line.startswith('#'):
            line = line.strip()
            dir_list.append(line)
    return dir_list


