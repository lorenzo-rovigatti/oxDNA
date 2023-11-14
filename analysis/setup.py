from setuptools import setup, find_packages
import re
import os
import numpy as np
from glob import glob

### Cython setup
USING_CYTHON = True
try:
    from Cython.Distutils.extension import Extension
    from Cython.Build import cythonize
except ImportError:
    from setuptools import Extension
    USING_CYTHON = False

OAT_DIR = os.path.dirname(__file__)

ext = 'pyx' if USING_CYTHON else 'c'
sources = glob('src/oxDNA_analysis_tools/cython_utils/*.{}'.format(ext))
extensions = [
    Extension('oxDNA_analysis_tools.UTILS.'+source.split(os.path.sep)[-1].split('.')[0],
              sources=[source]) for source in sources]

if USING_CYTHON:
    extensions = cythonize(extensions, compiler_directives={'language_level' : "3"})

# Get version number
def get_property(prop):
    init_location = os.path.join(OAT_DIR, 'src/oxDNA_analysis_tools/__init__.py')
    result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop), open(init_location).read())
    return result.group(1)

#invoke actual setup
setup(name="oxDNA-analysis-tools", version=get_property('__version__'), ext_modules=extensions, include_dirs=[np.get_include()])

#Notification about command line interface and autocompletes
print("\n\n################################################################################")
print("\nAll set up! All scripts in the package can now be run from the command line with\n\n    oat <script_name>\n\nimporting functions in your own scripts can be done with\n\n    from oxDNA_analysis_tools.<file name> import <function>")

print("")

print("If you would like to enable autocomplete in Bash, either add\n\n    source oat-completion.sh\n\nto your .bashrc or run \n\n    cat oat-completion.sh >> ~/.bash_completion\n\nto make autocompletes available to all users.")
