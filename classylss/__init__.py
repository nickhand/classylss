from .version import version as __version__

from .gcl import Cosmology
from .gcl import ClassEngine
from .gcl import ClassParams 
from .gcl import transfers
from .gcl import cltypes

from . import power
from . import tools
from . import correlation

def _init():
    """
    Set up the path of data files, which are installed to the package directory.
    """

    import os

    path = os.path.dirname(__file__)
    path = os.path.join(path, 'data')

    # setting static variables with swig is tricky.
    # see http://www.swig.org/Doc3.0/SWIGDocumentation.html#Python_nn20

    from .gcl import cvar

    cvar.ClassEngine_Alpha_inf_hyrec_file = os.path.join(path, 'hyrec', 'Alpha_inf.dat')
    cvar.ClassEngine_R_inf_hyrec_file = os.path.join(path, 'hyrec', 'R_inf.dat')
    cvar.ClassEngine_two_photon_tables_hyrec_file = os.path.join(path, 'hyrec', 'two_photon_tables.dat')
    cvar.ClassEngine_sBBN_file = os.path.join(path, 'bbn', 'sBBN.dat')

_init(); del _init