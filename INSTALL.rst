
*****************************
Installing the Ptools library
*****************************

.. contents::
    :backlinks: none

Quick setup guide
=================

This is the very minimal set of instructions required to install Ptools
in a Python virtual environment.
It assumes python 2.7 and all Ptools dependencies have been duly installed::

    $ virtualenv ptools-env
    $ source ptools-env/bin/activate
    (ptools-env) $ pip install cython
    (ptools-env) $ git clone https://github.com/ptools/ptools.git    
    (ptools-env) $ cd ptools
    (ptools-env) $ python setup.py install


For more details, see the full documented procedure described below.


Building and installing Ptools
==============================

The Ptools library has few dependencies namely :

- python 2.7.xx
- a C++ compiler (e.g. the `GNU Compiler`_),
- a developers version of the Python_ interpreter,
- the Boost_ library, 
- the f2c_ library,
- the Cython_ Python package,
- the Git_ version control system,
- the pip_ Python package will be required to install more Python packages.


Installing Ptools dependencies
------------------------------

GCC, Python & Git
^^^^^^^^^^^^^^^^^

Here is the set of instructions required to install all Ptools dependencies 
on a Debian-based system (including Ubuntu)::

    $ apt-get update
    $ apt-get install g++ python-dev python-pip git 

Those instructions can be easily adapted to other systems (e.g. RedHat)


Boost & F2C libraries
^^^^^^^^^^^^^^^^^^^^^

The standard way to install those libraries would be to use the system packages
in the same fashion as we did for ``gcc``, ``python`` and ``git``::

    $ apt-get install libboost-dev libf2c2-dev

Alternatively, you can use the Ptools approved legacy versions of 
the Boost and f2c libraries (see details in `Building Ptools with legacy Boost and f2c libraries`_.


Ptools itself
-------------

We recommand using a virtual environment. Virtual environments are usefull
to isolate Python packages from the rest of the Python packages installed
on your system. Hence, it dramatically limits the potential conflicts.
If ``virtualenv`` is not already present on your system, install it using ``pip``::

    $ pip install virtualenv

Assuming all dependencies are installed in standard locations, here is
the set of instructions required to build Ptools.

1. First you need to setup a virtual environment in which Ptools will be installed::

    $ virtualenv ptools-env
    $ source ptools-env/bin/activate

2. Install ``Cython``::

    (ptools-env) $ pip install cython

   Note that, as you are running in a virtual environment, from now
   Python packages that you install are not installed system-wide but
   only in your local environment.

3. Retrieve Ptools sources from its GitHub repository::

    (ptools-env) $ git clone https://github.com/ptools/ptools.git

4. Build and install Ptools::

    (ptools-env) $ cd ptools
    (ptools-env) $ python setup.py install


Using Ptools once it has been installed
=======================================

If you used ``virtualenv``, Ptools has been installed to your ``ptools-env``
environment. To use it in a new terminal, you need to activate the virtual
environment::

    $ source ptools-env/bin/activate

To check everything worked fine::

    (ptools-env) $ python -c 'import ptools'




Building Ptools with legacy Boost and f2c libraries
===================================================

First retrieve ptools sources (see `quick setup guide`_).
Then the procedure has been made quite straight forward by Ptools developers team::

    (ptools-env) $ cd ptools
    (ptools-env) $ python setup.py build_ext --use-legacy-f2c
    (ptools-env) $ python setup.py install


Troubleshooting
===============

boost/shared_array.hpp: No such file or directory
-------------------------------------------------

Ptools did not find the boost headers. You can explicitely tell Ptools where
to find it using the ``--with-boost-include-dir`` option at build time.
In this example, ``Boost`` has been installed in the directory ``/opt/boost``::

    (ptools-env) $ python setup.py build_ext --with-boost-include-dir=/opt/boost/include
    (ptools-env) $ python setup.py install

Alternatively, you can use the ``BOOST_INCLUDE_DIR`` environment variable::

    (ptools-env) $ export BOOST_INCLUDE_DIR=/opt/boost/include
    (ptools-env) $ python setup.py install


f2c.h: No such file or directory
--------------------------------

Ptools did not find the f2c headers. You can explicitely tell Ptools where
to find it using the ``--with-f2c-include-dir`` option at build time.
Importantly, this option is paired with the ``--with-f2c-library`` which
informs the location of the ``libf2c.a`` library.
In this example, ``f2c`` has been installed in the directory ``/opt/f2c``::

    (ptools-env) $ python setup.py build_ext --with-f2c-include-dir=/opt/f2c/include/ --with-f2c-library=/opt/f2c/lib/libf2c.a
    (ptools-env) $ python setup.py install

Alternatively, you can use the ``F2C_INCLUDE_DIR`` and ``F2C_LIBRARY``
environment variable::

    (ptools-env) $ export F2C_INCLUDE_DIR=/opt/f2c/include
    (ptools-env) $ export F2C_LIBRARY=/opt/f2c/lib/libf2c.a
    (ptools-env) $ python setup.py install


ImportError: [...]_ptools.so: undefined symbol: etime\_
-------------------------------------------------------

This error message occurs when importing ``ptools``. It is due to an error with
the linkage with the f2c library. To solve this problem, specify ``libf2c.a``
location as described in `f2c.h: No such file or directory`_.






.. _Boost: http://www.boost.org/
.. _f2c: http://www.netlib.org/f2c/
.. _Cython: http://cython.org/
.. _GNU Compiler: http://gcc.gnu.org/
.. _Python: http://www.python.org/
.. _Git: http://git-scm.com/
.. _pip: https://pypi.python.org/pypi/pip
