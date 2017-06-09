from .version import version as __version__

from .gcl import Cosmology
from .gcl import ClassEngine
from .gcl import ClassParams 
from .gcl import transfers
from .gcl import cltypes

from . import power
from . import tools
from . import correlation

def get_data_files():
    """
    Returns the path of data files, which are installed to the package directory.
    """

    import os

    path = os.path.dirname(__file__)
    path = os.path.join(path, 'data')
    r = dict(
        Alpha_inf_hyrec_file = os.path.join(path, 'hyrec', 'Alpha_inf.dat'),
        R_inf_hyrec_file = os.path.join(path, 'hyrec', 'R_inf.dat'),
        two_photon_tables_hyrec_file = os.path.join(path, 'hyrec', 'two_photon_tables.dat'),
        sBBN_file = os.path.join(path, 'bbn', 'sBBN.dat'),
    ) 
    return r

def _init():
    r = get_data_files()

    # setting static variables with swig is tricky.
    # see http://www.swig.org/Doc3.0/SWIGDocumentation.html#Python_nn20

    from .gcl import cvar

    cvar.ClassEngine_Alpha_inf_hyrec_file = r['Alpha_inf_hyrec_file']
    cvar.ClassEngine_R_inf_hyrec_file = r['R_inf_hyrec_file']
    cvar.ClassEngine_two_photon_tables_hyrec_file = r['two_photon_tables_hyrec_file']
    cvar.ClassEngine_sBBN_file = r['sBBN_file']

_init(); del _init
