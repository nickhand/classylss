from .version import version as __version__

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

def load_precision(precision):
    import os

    path = os.path.dirname(__file__)
    path = os.path.join(path, 'data', precision)

    r = dict()
    with open(path, 'r') as f:
        exec(f.read(), {}, r)

    return r

def load_ini(filename):
    """
    Read a CLASS ``.ini`` file, returning a dictionary of parameters

    Parameters
    ----------
    filename : str
        the name of the file to read
    """
    pars = {}

    with open(filename, 'r') as ff:

        # loop over lines
        for lineno, line in enumerate(ff):

            # skip any commented lines with #
            if '#' in line: line = line[line.index('#'):]

            # must have an equals sign to be valid
            if "=" not in line: continue

            # extract key and value pairs
            fields = line.split("=")
            assert len(fields) == 2, "error reading line number %d" %lineno
            pars[fields[0].strip()] = fields[1].strip()

    return pars
