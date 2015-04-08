import glob
import imp
import os
import sys

"""
Goal is to dynamically load classes in Test_*.py into UserTest
This is a horrible hack, but I've tried loads and it's better than eval works.

Alternatives...
    imp.find... and imp.load... don't put the module in the global namespace.
    eval("import %s" % test_name)
"""
# This hack is already bad enough, don't want folks importing *
__all__ = []

FILE_DIR = os.path.dirname(os.path.realpath(__file__))
TESTFILE_LIST = glob.glob("UserTest/Test*.py")
TEST_DICT = {}
for test_name in TESTFILE_LIST:
    test_name = test_name.strip('.py')
    # Attempt to append to globals/sys.modules fails
    name = os.path.basename(test_name)
   
    f, pathname, description = imp.find_module(test_name, [os.path.join(FILE_DIR, "..")])
    module = imp.load_module(name, f, pathname, description)
    test_class = module.__getattribute__(name)
    globals()[name] = module
    sys.modules[name] = module
    TEST_DICT.update({name: test_class})

    # This method works in Python2.4, but not Python2.7
    # Apparently importing a filepath this way was a bug
    # See: http://xkcd.com/1172/
    #module  = __import__(test_name).__getattribute__(name)
    #globals()[name] = module
    #sys.modules[name] = module

# Erase all trace of this horrible hack
del glob, imp, os, sys, description, f, pathname, module, name, test_class, test_name, TESTFILE_LIST

