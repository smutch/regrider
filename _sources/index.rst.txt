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

The easiest way to install and work with the code is using `Spack`_. Once you
have `Spack`_ installed::

    spack env create -d . spack.yaml
    spack env activate .

You can then build the code with the usual CMake commands::

    cmake -S. -Bbuild
    cmake --build build

If you don't want to use Spack then the following libraries must be installed
and the appropropriate flags passed to CMake so it can find them (typically
`CMAKE_PREFIX_PATH=...`).

* `HDF5`_ (with c++ and high-level libraries enabled)
* `FFTW3`_ (with openmp support)
* `fmt`_ (compiled with c++11 standard)
* `Criterion`_

.. _Spack: https://spack.readthedocs.io
.. _HDF5: https://www.hdfgroup.org/solutions/hdf5/
.. _FFTW3: http://www.fftw.org/
.. _fmt: https://fmt.dev/latest/index.html
.. _Criterion: https://criterion.readthedocs.io/en/master/

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
