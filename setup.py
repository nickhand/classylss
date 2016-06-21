from numpy.distutils.core import setup, Extension
from numpy.distutils.command.build_clib import build_clib
from numpy.distutils.command.build_ext import build_ext
from numpy.distutils.command.build import build
from numpy.distutils.command.install import install

from glob import glob
import os
import numpy
import shutil

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
MINOR = 0
MICRO = 1
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

DISTNAME = 'classy_lss'
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
    filename = os.path.join(os.path.dirname(__file__), 'classy_lss', 'version.py')
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
    args = (package_basedir, CLASS_VERSION, prefix)
    command = 'sh %s/depends/install_class.sh %s %s' %args
    
    ret = os.system(command)
    if ret != 0:
        raise ValueError("could not build CLASS v%s" %CLASS_VERSION)
    
class build_external_clib(build_clib):
    
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

    # Calls the default run command, then deletes the build area
    # (equivalent to "setup clean --all").
    def run(self):
        install.run(self)
        shutil.rmtree("build")
  

# setup FFTW for GCL library
fftw_info = {'include_dirs':[], 'library_dirs':[]}
if 'FFTW_INC' in os.environ:
    fftw_info['include_dirs'] = [os.environ['FFTW_INC']]
if 'FFTW_LIB' in os.environ:
    fftw_info['library_dirs'] = [os.environ['FFTW_LIB']]

# GCL extension 
gcl_info = {}
gcl_info['sources'] = list(glob("classy_lss/_gcl/cpp/*cpp"))
gcl_info['include_dirs'] = ['classy_lss/_gcl/include'] + fftw_info['include_dirs']
gcl_info['library_dirs'] = fftw_info['library_dirs']
libgcl = ('gcl', gcl_info)

    
sources = list(glob("classy_lss/_gcl/python/*.i")) + ['classy_lss/gcl.i']    
ext = Extension(name='classy_lss._gcl',
                sources=['classy_lss/gcl.i'],
                swig_opts=['-c++', '-Wall'], 
                extra_link_args=["-g", '-fPIC'],
                extra_compile_args=["-fopenmp", "-O2", '-std=c++11'],
                libraries=['gcl', 'class', 'gomp', 'fftw3'],
                runtime_library_dirs=fftw_info['library_dirs']
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
      py_modules = ["classy_lss.gcl"]
)

