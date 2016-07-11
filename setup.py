from numpy.distutils.core import setup, Extension
from numpy.distutils.command.build_clib import build_clib
from numpy.distutils.command.build_ext import build_ext
from numpy.distutils.command.build import build
from numpy.distutils.command.install import install

from glob import glob
import os
import numpy
import shutil

global install_dir

def get_latest_class_version():
    """
    Parse the ``class_public`` github API to determine the 
    latest version tag
    
    Returns
    -------
    name : str
        the latest version number
    """
    from distutils.version import LooseVersion
    import requests
    import json
    
    url = "https://api.github.com/repos/lesgourg/class_public/tags"
    req = requests.get(url)
    if req.status_code == requests.codes.ok:
        j = json.loads(req.text.encode(req.encoding))
        
        latest_version = LooseVersion("0.0.0")
        for release in j:
            version = release['name']
            if version.startswith("v"): version = version[1:]
            if LooseVersion(version) > latest_version:
                latest_version = LooseVersion(version)
        if latest_version == '0.0.0':
            raise ValueError("cannot find latest CLASS version to download")
        
        return str(latest_version)
    
    else:
        raise ValueError("cannot find latest CLASS version to download")

package_basedir = os.path.abspath(os.path.dirname(__file__))
CLASS_VERSION = get_latest_class_version()

MAJOR = 0
MINOR = 1
MICRO = 0
ISRELEASED = True
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

DISTNAME = 'classylss'
AUTHOR = 'Nick Hand'
AUTHOR_EMAIL = 'nicholas.adam.hand@gmail.com'
INSTALL_REQUIRES = ['numpy', 'astropy', 'requests']
DESCRIPTION = "python binding of CLASS for large-scale structure calculations"

if not ISRELEASED:
    VERSION += '.dev'
    
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

def build_class(prefix):
    """
    Download and build CLASS
    """
    # latest class version and download link       
    args = (package_basedir, CLASS_VERSION, prefix, install_dir)
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

        build_class(self.class_build_dir)
        link_objects = ['libclass.a']
        link_objects = list(glob(os.path.join(self.class_build_dir, '*', 'libclass.a')))
        
        self.compiler.set_link_objects(link_objects)
        self.compiler.library_dirs.insert(0, os.path.join(self.class_build_dir, 'lib'))        
        
        for (library, build_info) in libraries:
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
            
        build_ext.run(self)
        
        
class custom_install(install):
    """
    Custom install that deletes the ``build`` directory if successfull
    """
    def run(self):
        
        # set the global install dir
        global install_dir
        install_dir = os.path.join(self.install_lib, self.config_vars['dist_name'], 'class')
        install.run(self)
        shutil.rmtree("build")
  

gcl_sources = list(glob("classylss/_gcl/cpp/*cpp"))
fftlog_sources = list(glob("classylss/_gcl/extern/fftlog/*f"))

# GCL extension 
gcl_info = {}
gcl_info['sources'] =  gcl_sources + fftlog_sources 
gcl_info['include_dirs'] = ['classylss/_gcl/include']
gcl_info['language'] = 'c++'
libgcl = ('gcl', gcl_info)
    
sources = list(glob("classylss/_gcl/python/*.i")) + ['classylss/gcl.i']    
ext = Extension(name='classylss._gcl',
                sources=['classylss/gcl.i'],
                swig_opts=['-c++', '-Wall'], 
                extra_link_args=["-g", '-fPIC'],
                extra_compile_args=["-fopenmp", "-O2", '-std=c++11'],
                libraries=['gcl', 'class', 'gomp', 'gfortran']
                )


setup(name=DISTNAME,
      version=VERSION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      description=DESCRIPTION,
      install_requires=INSTALL_REQUIRES,
      ext_modules = [ext],
      libraries=[libgcl],
      cmdclass = {
          'build_clib': build_external_clib,
          'build_ext': custom_build_ext,
          'install': custom_install
      },
      py_modules = ["classylss.gcl"]
)

