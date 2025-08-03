Installation
============

MolR can be installed using pip or from source for development.

Requirements
------------

MolR requires Python 3.8 or later and the following dependencies:

- NumPy ≥ 1.20.0
- SciPy ≥ 1.7.0 (for KDTree spatial indexing)
- pyparsing ≥ 3.0.0 (for selection language)
- pdbreader ≥ 0.1.0 (for PDB file parsing)
- mmcif ≥ 0.1.0 (for mmCIF file parsing)

Installing from PyPI
--------------------

The simplest way to install MolR is using pip:

.. code-block:: bash

   pip install molr

This will install MolR and all its dependencies.

Installing from Source
----------------------

For development or to get the latest features, you can install from source:

.. code-block:: bash

   git clone https://github.com/abhishektiwari/molr.git
   cd molr
   pip install -e .

Development Installation
------------------------

If you plan to contribute to MolR or run tests, install with development dependencies:

.. code-block:: bash

   git clone https://github.com/abhishektiwari/molr.git
   cd molr
   pip install -r requirements-dev.txt

This installs additional tools for development:

- pytest (for running tests)
- black (for code formatting)
- mypy (for type checking)
- flake8 (for linting)
- isort (for import sorting)
- pytest-cov (for test coverage)

Virtual Environment Setup
-------------------------

We recommend using a Python virtual environment for development:

Using pyenv (recommended):

.. code-block:: bash

   pyenv virtualenv 3.10.0 molr
   pyenv activate molr
   pip install -r requirements-dev.txt

Using venv:

.. code-block:: bash

   python -m venv molr-env
   source molr-env/bin/activate  # On Windows: molr-env\\Scripts\\activate
   pip install -r requirements-dev.txt

Verifying Installation
----------------------

To verify your installation, you can run:

.. code-block:: python

   import molr
   print(molr.__version__)
   
   # Test basic functionality
   structure = molr.Structure(n_atoms=10)
   print(f"Created structure with {structure.n_atoms} atoms")

Running Tests
-------------

After development installation, you can run the test suite:

.. code-block:: bash

   # Run all tests except slow ones
   make test
   
   # Run all tests including slow ones
   make test-all
   
   # Run with coverage report
   make test-coverage

Building Documentation
----------------------

To build the documentation locally:

.. code-block:: bash

   cd docs
   make html

The documentation will be available at ``docs/build/html/index.html``.

Troubleshooting
---------------

**Import Error**: If you get import errors, ensure all dependencies are installed:

.. code-block:: bash

   pip install -r requirements.txt

**NumPy Version**: MolR requires NumPy 1.20.0 or later. Check your version:

.. code-block:: python

   import numpy
   print(numpy.__version__)

**SciPy KDTree**: If spatial queries are slow, ensure SciPy is properly installed:

.. code-block:: python

   from scipy.spatial import cKDTree
   print("SciPy KDTree available")

For more help, please open an issue on the `GitHub repository <https://github.com/abhishektiwari/molr>`_.