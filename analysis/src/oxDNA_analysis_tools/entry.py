# This script sets up the command line entry points to create the namespace "oat"

# Created by Erik Poppleton, October 2021

def main():
    from sys import argv
    from importlib import import_module

    try:
        script = argv.pop(1)
    except:
        from os import listdir, path
        print("Please select an analysis to run.  The available options are:")
        contents =listdir(path.dirname(__file__))

        #list the scripts that aren't __init__.py and entry.py
        contents = [c.replace('.py', '') for c in contents if ('.py' in c and not '__' in c and not path.basename(__file__) in c)]
        contents.sort()
        print('\n'.join(contents))
        exit(1)

    to_call = import_module('oxDNA_analysis_tools.{}'.format(script))

    to_call.main()