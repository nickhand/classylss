from distutils.core import Command
from numpy.distutils.core import Extension
from numpy.distutils.command.build_clib import build_clib
from numpy.distutils.command.build_ext import build_ext
from numpy.distutils.command.build import build
from distutils.command.clean import clean

from glob import glob
import os
import numpy
import shutil

# use GNU compilers by default  
os.environ.setdefault("CXX", "g++")
os.environ.setdefault("F90", "gfortran")

# base directory of package
package_basedir = os.path.abspath(os.path.dirname(__file__))

# the CLASS version to install
CLASS_VERSION = '2.5.0'

MAJOR = 0
MINOR = 1
MICRO = 8
ISRELEASED = True
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

DISTNAME = 'classylss'
AUTHOR = 'Nick Hand'
AUTHOR_EMAIL = 'nicholas.adam.hand@gmail.com'
INSTALL_REQUIRES = ['numpy', 'astropy']
DESCRIPTION = "python binding of CLASS for large-scale structure calculations"
URL = "http://github.com/nickhand/classylss"

if not ISRELEASED: VERSION += '.dev0'

def write_version_py():
    cnt = """\
version = '%s'
class_version = '%s'
"""
    filename = os.path.join(os.path.dirname(__file__), 'classylss', 'version.py')
    a = open(filename, 'w')
    try:
        a.write(cnt % (VERSION, CLASS_VERSION))
    finally:
        a.close()
        
write_version_py()
    
def check_swig_version():
    """
    Check the version of swig, >= 3.0 is required
    """
    import subprocess, re
    try:
        output = subprocess.check_output(["swig", "-version"])
    except OSError:
        raise ValueError("`swig` not found on PATH -- installation cannot proceed")
        
    try:
        version = re.findall("SWIG Version [0-9].[0-9].[0-9]", output)[0].split()[-1]
    except:
        return
        
    # need >= 3.0
    if version < "3.0":
        raise ValueError("the version of `swig` on PATH must greater or equal to 3.0")
    
def build_CLASS(prefix):
    """
    Download and build CLASS
    """
    # latest class version and download link       
    args = (package_basedir, CLASS_VERSION, prefix, "/opt/class/willfail")
    command = 'sh %s/depends/install_class.sh %s %s %s' %args
    
    ret = os.system(command)
    if ret != 0:
        raise ValueError("could not build CLASS v%s" %CLASS_VERSION)
    
class build_external_clib(build_clib):
    """
    Custom command to build CLASS first, and then GCL library
    """
    def finalize_options(self):
        
        build_clib.finalize_options(self)    
        
        # create the CLASS build directory and save the include path
        self.class_build_dir = self.build_temp
        self.include_dirs.insert(0, os.path.join(self.class_build_dir, 'include'))

    def build_libraries(self, libraries):

        # build CLASS first
        build_CLASS(self.class_build_dir)

        # update the link objects with CLASS library
        link_objects = ['libclass.a']
        link_objects = list(glob(os.path.join(self.class_build_dir, '*', 'libclass.a')))
                        
        self.compiler.set_link_objects(link_objects)
        self.compiler.library_dirs.insert(0, os.path.join(self.class_build_dir, 'lib'))        
        
        for (library, build_info) in libraries:
            
            # check swig version
            if library == "gcl": check_swig_version()
            
            # update include dirs
            self.include_dirs += build_info.get('include_dirs', [])
        
        build_clib.build_libraries(self, libraries)

class custom_build_ext(build_ext):
    """
    Custom extension building to grab include directories
    from the ``build_clib`` command
    """
    def finalize_options(self):
        build_ext.finalize_options(self)
        self.include_dirs.append(numpy.get_include())

    def run(self):
        if self.distribution.has_c_libraries():
            self.run_command('build_clib')
            build_clib = self.get_finalized_command('build_clib')
            self.include_dirs += build_clib.include_dirs
            self.library_dirs += build_clib.compiler.library_dirs
            
        # copy data files from temp to classlssy package directory
        shutil.rmtree(os.path.join(self.build_lib, 'classylss', 'data'), ignore_errors=True)
        shutil.copytree(os.path.join(self.build_temp, 'data'), os.path.join(self.build_lib, 'classylss', 'data'))
            
        build_ext.run(self)
        
class custom_clean(clean):

    def run(self):

        # run the built-in clean
        clean.run(self)
            
        # remove the CLASS tmp directories
        os.system("rm -rf depends/tmp*")
        
gcl_sources = list(glob("classylss/_gcl/cpp/*cpp"))
fftlog_sources = list(glob("classylss/_gcl/extern/fftlog/*f"))

# GCL extension 
gcl_info = {}
gcl_info['sources'] =  gcl_sources + fftlog_sources 
gcl_info['include_dirs'] = ['classylss/_gcl/include']
gcl_info['language'] = 'c++'
gcl_info['extra_compiler_args'] = ["-fopenmp", "-O2", '-std=c++11']
libgcl = ('gcl', gcl_info)
    
sources = list(glob("classylss/_gcl/python/*.i")) + ['classylss/gcl.i']    
ext = Extension(name='classylss._gcl',
                sources=['classylss/gcl.i'],
                depends=['classylss/_gcl/python/*.i'],
                swig_opts=['-c++', '-Wall'], 
                extra_link_args=["-g", '-fPIC'],
                extra_compile_args=['-fopenmp'],
                libraries=['gcl', 'class', 'gomp', 'gfortran']
                )
    
if __name__ == '__main__':
    
    from numpy.distutils.core import setup
    setup(name=DISTNAME,
          version=VERSION,
          author=AUTHOR,
          author_email=AUTHOR_EMAIL,
          description=DESCRIPTION,
          url=URL,
          requires=INSTALL_REQUIRES,
          ext_modules = [ext],
          libraries=[libgcl],
          cmdclass = {
              'build_clib': build_external_clib,
              'build_ext': custom_build_ext,
              'clean': custom_clean,
          },
          py_modules = ["classylss.gcl"]
    )
