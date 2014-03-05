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

class Test_MeanPE(BaseTest.Tester):
    """
    Instance of Tester class for testing the Mean PE
    Pass if meanPE (output_files) +/- 5*SEM (expected_result)
    """
    def do_test(self):
        """
        Extend do_test
        """
        energy_filepath = self.obj_runner.run_input_options['energy_file']
        if BaseTest.check_file(energy_filepath, sys_exit=False):
            energy = np.loadtxt(energy_filepath)
        else:
            self.test_pass = False
            self.test_pass_info = "Test Failed. File not found: " + energy_filepath
            return self.test_pass

        # PE: Potential Energy
        PE = energy[:,1]
        meanPE = np.mean(PE)
        stdPE = np.std(PE)
        semPE = stdPE/np.sqrt(np.size(PE))
        self.test_results_dict.update({'meanPE':meanPE, 'stdPE':stdPE, 'semPE':semPE})

        if not self.expected_results_dict:
            self.test_pass = False
            self.test_pass_info = "Test Failed. Empty comparison file: " + self.expected_results_filepath
            return self.test_pass

        if meanPE < float(self.expected_results_dict['meanPE']) + 5*float(self.expected_results_dict['semPE']) and \
                meanPE > float(self.expected_results_dict['meanPE']) - 5*float(self.expected_results_dict['semPE']):
                    self.test_pass = True
                    self.test_pass_info = "Test Passed"
        return self.test_pass
