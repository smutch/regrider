.. Regrider documentation master file, created by
   sphinx-quickstart on Fri Dec 20 00:10:04 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Regrider's documentation!
=========================

Downsample gbpTrees and VELOCIraptor trees using FFTW.

.. .. toctree::
..    :maxdepth: 2
..    :caption: Contents:

Usage
=====

.. code-block:: man

   Usage:
     regrider [OPTION...]
   
     -d, --dim arg           new grid dimension
     -g, --gbptrees arg      input gbpTrees grid file
     -v, --velociraptor arg  input VELOCIraptor grid file
     -o, --output arg        output file name
     -h, --help              show help


API
===

.. doxygenfile:: gbptrees.hpp

.. doxygenclass:: Grid
   :members:


Indices and tables
==================

* :ref:`genindex`
