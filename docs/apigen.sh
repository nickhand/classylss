# bash

if ! python -c 'import numpydoc'; then pip install --user numpydoc; fi
if ! python -c 'import sphinx'; then pip install --user sphinx; fi

pushd ..
python setup.py build_ext --inplace
popd

sphinx-apidoc -H "API Reference" -e -f -o . ../classylss ../classylss/tests
