classy_lss
----------
python binding of CLASS for large-scale structure calculations

Dependencies
------------

The latest version of CLASS will automatically be downloaded and linked against. FFW3 is required, and the
user should specify its location as detailed in the next section.


Installation
------------

The packaged can be cloned from github as

.. code:: bash

    git clone https://github.com/nickhand/classy_lss.git
    cd classy_lss
    
If the ``fftw3`` is installed in a non-standard location, the user should specify the
location using environment variables:

.. code:: bash

    export FFT_INC=/path/to/fftw/include
    export FFT_DIR=/path/to/fftw/lib

and then the package can be installed using

.. code:: bash
    
    export CFLAGS='-std=c++11 -fopenmp'; pip install .
    
This procedue has been tested on Mac and Linux machines.
    
