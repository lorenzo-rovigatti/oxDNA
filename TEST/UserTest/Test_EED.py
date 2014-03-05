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

class Test_EED(BaseTest.Tester):
    """
    Instance of Tester class for testing the End-to-End-Distance.
    Pass if meanEED (output_files) +/- 5*SEM (expected_result)
    """
    def do_test(self):
        """
        Extend do_test from class Tester
        """
        run_trajectory = self.obj_runner.run_trajectory
        strand_length = run_trajectory[0]._strands[0].get_length()
        mean_dist_theory = np.sqrt(strand_length)
        mean_dist_error = 0.1 * mean_dist_theory

        dist_list = []
        for traj in run_trajectory:
            nt_first = traj._nucleotides[0]
            nt_last  = traj._nucleotides[-1]
            dr = nt_first.distance(nt_last, box=traj._box)
            dist = np.sqrt(np.dot(dr,dr))
            dist_list.append(dist)
        dist_list = np.array(dist_list)
        mean_dist = np.mean(dist_list)
        std_dist = np.std(dist_list)
        sem_dist = std_dist/np.sqrt(np.size(dist_list))

        self.test_results_dict.update({'mean_dist':mean_dist, 'std_dist':std_dist, 'sem_dist':sem_dist})

        if not self.expected_results_dict:
            self.test_pass = False
            self.test_pass_info = "Test Failed. Empty comparison file: " + self.expected_results_filepath
            return self.test_pass

        if mean_dist < float(self.expected_results_dict['mean_dist']) + 5*float(self.expected_results_dict['sem_dist']) and \
                mean_dist > float(self.expected_results_dict['mean_dist']) - 5*float(self.expected_results_dict['sem_dist']):
                    self.test_pass = True
                    self.test_pass_info = "Test Passed"
        return self.test_pass
