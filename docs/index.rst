.. Regrider documentation master file, created by
   sphinx-quickstart on Fri Dec 20 00:10:04 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Regrider
========

| **Downsample gbpTrees and VELOCIraptor trees using FFTW.**
| https://github.com/smutch/regrider


Installation
============

Coming soon...


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

A utility script is also provided to downsample a directory of VELOCIraptor grids:

.. code-block:: man

    Usage:
        ./downsample_vr_directory.sh [OPTION...]

        -b <dir>    Bin directory of hdf5 (location of h5ls and h5copy)
        -s <dir>    Source directory containing VR grids
        -d <dir>    Destination directory to place VR grids (will be created if necessary)
        -n <int>    New grid dimension (n x n x n)

    All options must be present.


.. toctree::
   :maxdepth: 2
   :caption: Contents

   GBPtrees <gbptrees>
   VELOCIraptor <velociraptor>
   grid
   utils


Indices and tables
==================

* :ref:`genindex`
