Installation
============

From source (current)
---------------------

Clone the repository and install with pip::

    git clone https://github.com/dehaenw/lacan.git
    cd lacan
    pip install .

This installs LACAN in your current Python environment along with its
dependencies. 

Requirements
------------

* Python ≥ 3.9
* RDKit

PyPI
----

Once available, LACAN can also be installed directly from PyPI::

    pip install lacan

Development install
-------------------

To run the test suite or contribute to the code, install in editable mode with
the development extras::

    pip install -e ".[dev]"
    pytest          # run the full test suite (~155 tests)
