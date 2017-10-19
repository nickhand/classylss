from .version import version as __version__, class_version

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

def _find_file(filename):
    """
    Find the file path, first checking if it exists and then looking in the
    data directory
    """
    import os
    if os.path.exists(filename):
        path = filename
    else:
        path = os.path.dirname(__file__)
        path = os.path.join(path, 'data', filename)

    if not os.path.exists(path):
        raise ValueError("cannot locate file '%s'" %filename)

    return path

def load_precision(filename):
    """
    Load a CLASS precision file into a dictionary.

    Parameters
    ----------
    filename : str
        the name of an existing file to load, or one in the files included
        as part of the CLASS source

    Returns
    -------
    dict :
        the precision parameters loaded from file
    """
    # also look in data dir
    path = _find_file(filename)

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
        the name of an existing parameter file to load, or one included as
        part of the CLASS source

    Returns
    -------
    dict :
        the input parameters loaded from file
    """
    # also look in data dir
    path = _find_file(filename)

    pars = {}
    with open(path, 'r') as ff:

        # loop over lines
        for lineno, line in enumerate(ff):
            if not line: continue

            # skip any commented lines with #
            if '#' in line: line = line[line.index('#')+1:]

            # must have an equals sign to be valid
            if "=" not in line: continue

            # extract key and value pairs
            fields = line.split("=")
            if len(fields) != 2:
                import warnings
                warnings.warn("skipping line number %d: '%s'" %(lineno,line))
                continue
            pars[fields[0].strip()] = fields[1].strip()

    return pars
