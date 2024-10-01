import os
import re
import sys
import sysconfig
import platform
import subprocess
import multiprocessing

from distutils.version import LooseVersion
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.install_lib import install_lib

PACKAGE_NAME = "oxpy"


class CMakeExtension(Extension):

    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)
        
        
class InstallCMakeLibs(install_lib):
    """
    Get the libraries from the parent distribution, use those as the outfiles

    Skip building anything; everything is already built, forward libraries to
    the installation step
    """

    def run(self):
        """
        Copy libraries from the bin directory and place them as appropriate
        """

        self.announce("Moving library files", level=3)

        # We have already built the libraries in the previous build_ext step
        self.skip_build = True

        # Folder where the `oxpy` package has been placed by cmake. It is used by self.install
        self.build_dir = self.distribution.lib_dir
        self.outfiles = self.install()
        
        # I have copied this bit from the parent class
        if self.outfiles is not None:
            # always compile, in case we have any extension stubs to deal with
            self.byte_compile(self.outfiles)
            
    def get_outputs(self):
        """
        Overrides the parent class' method. Returns a list of the files copied over by the `run` method 
        """
        return self.outfiles


class CMakeBuild(build_ext):
    user_options = build_ext.user_options + [
        ('cuda', None, 'enable CUDA support')
    ]
    
    def initialize_options(self):
        build_ext.initialize_options(self)
        self.cuda = None

    def finalize_options(self):
        build_ext.finalize_options(self)

    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " + ", ".join(e.name for e in self.extensions))

        self.cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
        
        if platform.system() == "Windows":
            if self.cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)
            
    def build_extension(self, ext):
        self.announce("Preparing the build environment", level=3)
        
        cmake_args = []

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]
        native_generator_args = ['--']

        if platform.system() == "Windows":
            if sys.maxsize > 2 ** 32:
                cmake_args += ['-A', 'x64']
            native_generator_args += ['/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
                
        cmake_args += ['-DPython=On']
        if self.cuda is not None:
            cmake_args += ['-DCUDA=On']
        
        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            
            try:
                cpu_cores = int(os.getenv('SLURM_NTASKS'))
            except:
                cpu_cores = int(multiprocessing.cpu_count() / 2)
            
            if self.cmake_version < "3.14.0":
                native_generator_args += ["-j{}".format(cpu_cores)]
            else:
                build_args += ["-j {}".format(cpu_cores)]
                
        build_args += native_generator_args

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
            
        self.distribution.lib_dir = os.path.join(self.build_temp, "oxpy/python")
            
        self.announce("Configuring cmake project", level=3)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp)
        
        self.announce("Building the library", level=3)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

        self.announce("Compilation done", level=3)


setup(
    name = PACKAGE_NAME,
    use_scm_version = {
        "fallback_version": "3.3",
        },
    packages=find_packages(),
    setup_requires = ['setuptools-scm'],
    install_requires = [
        f"oxDNA_analysis_tools @ file://localhost/{os.getcwd()}/analysis/"
    ],
    dependency_links = [''],
    author = 'Lorenzo Rovigatti, Flavio Romano, Petr Sulc and others',
    author_email = 'lorenzo.rovigatti@uniroma1.it',
    description = 'A code primarily aimed at DNA and RNA coarse-grained simulations',
    long_description = open("./README.md", 'r').read(),
    long_description_content_type = "text/markdown",
    license = 'GNU GPL 3.0',
    ext_modules = [CMakeExtension('oxpy')],
    # add custom build_ext command
    cmdclass={
        'build_ext': CMakeBuild,
        'install_lib': InstallCMakeLibs,
        },
    zip_safe = False,
)
