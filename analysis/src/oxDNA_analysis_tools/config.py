#!/usr/bin/env python

from os import path
from typing import List
from sys import stderr
import argparse
from oxDNA_analysis_tools.UTILS.logger import log, logger_settings

#checking dependencies to make sure everything is correct
def check(to_check:List[str]=["python", "numpy", "matplotlib", "sklearn", "oxpy"]):
    """
        Checks if the dependencies are installed.

        Parameters:
            to_check (List[str]): list of package names to check
    """
    flag = False
    dependencies = {
        "numpy": 1.14,
        "matplotlib": 3.0,
        "sklearn": 0.21,
        "oxpy": 3.2,
    }
    real_names = {
        "numpy": "Numpy",
        "matplotlib": "MatPlotLib",
        "sklearn": "SciKit-Learn",
        "oxpy": "oxpy"
    }
    websites = {
        "numpy": "numpy.org", 
        "matplotlib": "matplotlib.org",
        "sklearn": "scikit-learn.org",
        "oxpy": "https://github.com/lorenzo-rovigatti/oxDNA"
    }

    #get version of this package
    oat = __import__('oxDNA_analysis_tools')
    log("oxDNA_analysis_tools version: {}".format(oat.__version__))
    log(f"running config.py installed at: {path.realpath(__file__)}")


    #check python version
    if "python" in to_check:
        from sys import version_info
        ver = '.'.join([str(i) for i in version_info[0:2]])
        log("Python version: {}".format('.'.join([str(i) for i in version_info[0:3]])))
        if version_info < (3, 9):
            flag = True
            log("Some scripts will not run with Python versions earlier than 3.9.  You have {}, please update your environment".format(version_info), level="warning")

    #check packages
    for package in to_check:
        if package == "python": continue
        try:
            mod = __import__(package)
            log("Package {} found. Version: {}".format(real_names[package], mod.__version__))
        except:
            flag = True
            # Log doesn't handle errors, but don't raise anything here because we want to continue execution
            print("ERROR: Unable to find package {}.  Please check your environment or follow the installation instructions at {}".format(real_names[package], websites[package]), file=stderr)
            continue
        ver = float('.'.join(mod.__version__.split(".")[0:2]))
        if ver < dependencies[package]:
            flag = True
            log("Your version for package {} is {}.  This tool was tested using {}.  You may need to update your environment".format(real_names[package], ver, dependencies[package]), level="warning")

    # Check for numpy header error
    try:
        from oxDNA_analysis_tools.UTILS.get_confs import cget_confs
    except ValueError as e:
        flag = True
        raise RuntimeError(f"Importing Cython file reading module failed with error:\n{e}\n\
                           This is generally because your version of Numpy is behind the version downloaded by pip during installation.\n\
                           Please try reinstalling oxDNA_analysis_tools with --no-build-isolation to use your local numpy header. Or update your numpy version.\n\
                           This can be done by navigating to the oxDNA/analysis directory and running:\n\
                           \t python -m pip install . --no-build-isolation")
    
    if flag:
        log("Some packages need to be installed/updated.", level="warning")
    else:
        log("No dependency issues found.")

    return flag

def set_chunk_size(chunk_size:int):
    """
        Sets the number of confs to read at a time for analyses.  This value is persistent between analyses.

        Parameters:
            chunk_size (int): number of confs to read at a time

        This will update a file called chunksize.py in the UTILS directory.
    """
    with open(path.realpath(__file__).strip('config.py')+"UTILS/chunksize.py", 'w') as f:
        f.write("CHUNKSIZE = "+str(chunk_size))

def get_chunk_size():
    """
        Gets the current chunk size.

        Returns:
            int: number of confs to read at a time
    """
    try:
        from oxDNA_analysis_tools.UTILS.chunksize import CHUNKSIZE
        log("Analyses will be computed in chunks of {} configurations at a time".format(CHUNKSIZE))
        log("You can modify this number by running oat config -n <number>, which will be persistent between analyses.")
    except:
        raise Exception("Unable to read chunksize from file. UTILS/chunksize.py should contain a line like CHUNKSIZE = 100")

def cli_parser(prog="config.py"):
    parser = argparse.ArgumentParser(prog = prog, description='Configure oxDNA_analysis_tools')
    parser.add_argument('-n', '--chunk_size', type=int, help='Number of configurations to per chunk.  Persistent across analyses.')
    parser.add_argument('-q', metavar='quiet', dest='quiet', action='store_const', const=True, default=False, help="Don't print 'INFO' messages to stderr")
    return parser

def main():
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()
    
    logger_settings.set_quiet(args.quiet)
    if args.chunk_size:
        set_chunk_size(args.chunk_size)
        log("Future analyses will calculate in blocks of {} confs at a time".format(args.chunk_size))

    check(["python", "numpy", "matplotlib", "sklearn", "oxpy"])

    print()
    get_chunk_size()

if __name__ == '__main__':
    main()
    
