import logging
import inspect
from optparse import OptionParser
import os
import sys

import BaseTest
import UserTest


if __name__ == '__main__':
    """
    Execute intergration testing for oxDNA using Runner, Tester, TestSuite framework implemented in BaseTest

    RunTest.py relies on three special files. These are not hardcoded into the BaseTest framework, which retains much more flexibility.
    INPUT --- oxDNA INPUT file for each Runner
    TEST_LIST --- File with a list of UserTest (specified by Tester module/class name).
    RUNNER_LIST --- File with a list of directories to search for INPUT and TEST_LIST files

    Order of Battle:
    1) Search through all subdirectories of directories in the RUNNER_LIST to find INPUT/TEST_LIST files.
    2) Create Runners for each INPUT file. Each TestSuite consist of only one Runner.
    3) Run UserTest, designated in TEST_LIST, for each INPUT file
    """
    usage = "usage: %prog [options] <RUNNER_LIST>"
    parser = OptionParser(usage=usage)
    parser.add_option("--input_filename", dest="input_filename", default="INPUT",
            help="oxDNA input filename", metavar="FILE")
    parser.add_option("--exec_filepath", dest="exec_filepath", default="../bin/oxDNA",
            help="oxDNA executable", metavar="FILE")
    parser.add_option("--base_path", action="store", dest="base_path",
            default=".", help="Set the base_path")
    parser.add_option("--test_list_filename", dest="test_list_filename", default="TEST_LIST",
            help="Select which Tests to run", metavar="FILE")
    parser.add_option("--testclass_path", action="store", dest="testclass_path",
            default="UserTest", help="Set the path where UserTestes are stored")
    #parser.add_option("-d", "--dot_avg", action="store_true", dest="bool_dot_avg",
    #        default=False, help="Output average dot product instead of kinks")
    #parser.add_option("-i", "--increment", action="store", dest="increment",
    #        default=2, help="Set integer dot(n, n+INCREMENT), must be >= 1")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
            default=False, help="Verbose output")
    (options, args) = parser.parse_args()

    if len(args) > 1:
        parser.error("Too many args")
    elif len(args) == 1:
        runner_list_filename = str(args[0])
    else:
        runner_list_filename = 'RUNNER_LIST_BUILD'

    # Clean OptionParser inputs
    input_filename = str(options.input_filename)
    exec_filepath = str(options.exec_filepath)
    test_list_filename = str(options.test_list_filename)
    if BaseTest.check_path(options.base_path):
        base_path = str(options.base_path)
    if BaseTest.check_path(options.testclass_path):
        testclass_path = str(options.testclass_path)

    # Setup logging
    verbose = bool(options.verbose)
    logger = logging.getLogger('RunTest')
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARNING)

    # Load RUNNER_LIST from directory base_path
    dir_list = BaseTest.load_file(runner_list_filename)
    logger.info("Searching for oxDNA INPUT files in subdirectories of %s:\n%s", runner_list_filename, dir_list)

    # Walk to find INPUT in subdir of directories listed in RUNNER_LIST
    input_filepath_list = []
    for dir_name in dir_list:
        for dirpath, dirnames, filenames in os.walk(dir_name):
            if input_filename in filenames:
                input_filepath = os.path.join(dirpath, input_filename)
                if BaseTest.check_file(input_filepath, sys_exit=True):
                    input_filepath_list.append(input_filepath)
                    logger.info("oxDNA INPUT file (%s)", input_filepath)
                else:
                    logger.info("Skipping input file (%s)", input_filepath)

    #Used UserTest/__init__.py to bypass the modules. If I didn't, this code would work.
    # Use only those designated in TEST_LIST (same directory as oxDNA INPUT file by default)
    #tester_module_dict = dict(inspect.getmembers(sys.modules['UserTest'], inspect.ismodule))
    #tester_class_dict = {}
    #for key,value in tester_module_dict.items():
    #    tester_class_dict.update({key: value.__getattribute__(key)})

    # Dictionary of all potential UserTest
    #tester_class_dict = dict(inspect.getmembers(sys.modules['UserTest'], inspect.isclass))
    tester_class_dict = UserTest.TEST_DICT
    logger.info("Tester Classes Loaded:\n%s", tester_class_dict.keys())

    N_failed = 0
    N_total = 0
    for input_filepath in input_filepath_list:
        input_path, input_filename = os.path.split(input_filepath)
        runner_instance_list = [BaseTest.Runner(exec_filepath=exec_filepath,input_filename=input_filename, input_path=input_path, verbose=verbose)]

        # Load TEST_LIST from same directory as INPUT
        test_list_filepath = os.path.join(input_path, test_list_filename)

        # Check the TEST_LIST
        if not BaseTest.check_file(test_list_filepath, sys_exit=True):
            logger.info("Skipping INPUT (%s). No TEST_LIST found at path (%s)", input_filepath, test_list_filepath)
            continue
        test_list = BaseTest.load_file(test_list_filepath)

        tester_class_list = []
        for test_name in test_list:
            try:
                tester_class_list.append(tester_class_dict[test_name])
            except (KeyError):
                logger.info("Tester not found (%s). Skipping test for Runner at path (%s)", test_name, input_filepath)
                continue

        # Add tests to new TestSuite
        TS = BaseTest.TestSuite(verbose=verbose)
        for tester_class in tester_class_list:
            for runner_instance in runner_instance_list:
                obj_test = tester_class(runner_instance,
                        output_files_path=os.path.join(input_path, 'output_files/'),
                        expected_results_path=os.path.join(input_path, 'expected_results/'),
                        verbose=False)
                TS.add_test(obj_test)
        TS.run()
        TS.result()
        N_failed += TS.N_failed
        N_total += TS.N_total

    # Print overall summary at end
    print "\nFAILED/TOTAL: %i/%i\n" % (N_failed, N_total)

