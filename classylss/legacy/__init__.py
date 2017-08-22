from ..gcl import Cosmology
from ..gcl import ClassEngine
from ..gcl import ClassParams
from ..gcl import transfers
from ..gcl import cltypes

def _init():
    from .. import get_data_files
    r = get_data_files()

    # setting static variables with swig is tricky.
    # see http://www.swig.org/Doc3.0/SWIGDocumentation.html#Python_nn20

    from ..gcl import cvar

    cvar.ClassEngine_Alpha_inf_hyrec_file = r['Alpha_inf_hyrec_file']
    cvar.ClassEngine_R_inf_hyrec_file = r['R_inf_hyrec_file']
    cvar.ClassEngine_two_photon_tables_hyrec_file = r['two_photon_tables_hyrec_file']
    cvar.ClassEngine_sBBN_file = r['sBBN_file']

_init(); del _init
