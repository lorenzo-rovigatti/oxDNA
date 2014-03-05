import logging
import sys

# Logging
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)

# Gracefully load numpy
try:
    import numpy as np
except (ImportError):
    logging.critical("Cannot import numpy. Skipping Test (%s)", __name__)

# Gracefully load BaseTest
try:
    sys.path.append('../')
    import BaseTest
except (ImportError):
    logging.critical("Cannot import BaseTest. Dying gracefully (%s)", __name__)
    sys.exit()

class Test_Run(BaseTest.Tester):
    """
    Instance of Tester class for testing if oxDNA ran.
    Quick and dirty!
    """
    def do_test(self):
        """
        Extend do_test from class Tester
        """
        if self.obj_runner.run_success:
            self.test_pass = True
            self.test_pass_info = "Test Passed"
        else:
            self.test_pass = False
            self.test_pass_info = "Test Failed. oxDNA Run Failed"
        self.test_results_dict.update({'run_success': self.test_pass})
        return self.test_pass
