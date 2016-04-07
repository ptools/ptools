# -*- coding: utf-8 -*-
import os
import subprocess
import sys
import textwrap

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


# !!! PLEASE override the following two variables
# !!! to define a custom path to required dependencies:

# user-defined path to libf2c.a
user_path_libf2c = ""

# user-defined path to f2c.h
user_path_f2c_h = ""

# user-defined path to boost library:
user_path_boost = ""


# =============================================================================

NOT_FOUND_MESSAGE = 'Note: You can add a custom search path by editing this '\
                    'file (setup.py).\nYou can also install locally the '\
                    'missing dependencies by running: sh ./install-deps.sh'


def git_version():
    """Return the git revision as a string."""
    cmd = ['git', 'show', '-s', '--format=%h %ci', 'HEAD']
    try:
        git_revision = subprocess.check_output(cmd).strip()
    except OSError:
        git_revision = 'Unknown version. Please use git to download PTools '\
                       'and get reliable versioning informations'
    return git_revision


def write_version_h(filename):
    """Write a header file with the current git revision hash."""
    git_revision = git_version()
    if git_revision.startswith('Unknown'):
        s = "WARNING: it seems that you don't have a git directory."\
            "While the library will compile correcly, informations about"\
            "the current ptools version will be missing. Please use git to"\
            "download PTools and get reliable versioning informations."
        print textwrap.fill(s)

    content = """
/*
 * This file was generated automatically.
 * You should not modify it manually, as it may be re-generated.
 */

#ifndef GITREV_H
#define GITREV_H
#define GIT_REVID   "%(git_revision)s"
#endif /* GITREV_H */
"""
    with open(filename, 'w') as f:
        f.write(content % {'git_revision': git_revision})


# == Methods to locate headers and libraries ==

def find_file(name, paths):
    """Try to locate a file in a given set of directories.

    Args:
        name(str): file to look for.
        paths(list[str]): directories to scan.

    Return:
        str: the first directory in which one file has been found
            or an empty string if no file has been found.
    """
    for p in paths:
        fullfilepath = os.path.join(p, name)
        if os.path.exists(fullfilepath):
            return fullfilepath
    return ''


def find_header(names, paths, useEnvPath=False):
    """Try to locate a file in a given set of directories.

    Args:
        names(list[str]): files to look for.
        paths(list[str]): directories to scan.

    Return:
        str: the first directory in which one file has been found
            or an empty string if no file has been found.
    """
    for p in paths:
        if os.path.exists(p):
            for n in names:
                if os.path.exists(os.path.join(p, n)):
                    return p
    return ''


def find_boost():
    """Try to locate the boost libraries (actually just look for
    the shared_array.hpp header file."""
    print "Trying to locate the boost libraries."
    boostdir = find_header(["boost/shared_array.hpp"],
                           [user_path_boost, "./ptools_dep/boost_1_55_0",
                            "/usr/include", "/opt/local/include"])
    if not boostdir:
        print "Cannot locate boost library", NOT_FOUND_MESSAGE
        sys.exit(1)
    return boostdir


def find_f2c():
    """Try to locate the f2c library and header."""
    print "Trying to locate the libf2c.a static library and f2c.h header."
    f2clib = find_file("libf2c.a",
                       [user_path_libf2c, "./ptools_dep/f2c", "/usr/lib/",
                        "/usr/local/lib/", "/usr/lib64/"])
    f2c_header = find_header(["f2c.h"],
                             [user_path_f2c_h, "./ptools_dep/f2c",
                              "/usr/include", "/usr/local/f2c/"])
    if not f2clib:
        print "Cannot locate libf2c", NOT_FOUND_MESSAGE
        sys.exit(1)
    return f2clib, f2c_header


def setup_package():
    boost_include_dir = find_boost()
    f2clib, f2c_include_dir = find_f2c()
    print "using boost from", boost_include_dir
    print "using f2clib from", f2clib

    write_version_h('headers/gitrev.h')

    sources = ['src/cython_wrappers.cpp',
               'src/atom.cpp',
               'src/attractrigidbody.cpp',
               'src/coordsarray.cpp',
               'src/mcopff.cpp',
               'src/rigidbody.cpp',
               'src/surface.cpp',
               'src/atomselection.cpp',
               'src/basetypes.cpp',
               'src/forcefield.cpp',
               'src/pairlist.cpp',
               'src/rmsd.cpp',
               'src/version.cpp',
               'src/attractforcefield.cpp',
               'src/coord3d.cpp',
               'src/geometry.cpp',
               'src/pdbio.cpp',
               'src/superpose.cpp',
               'src/scorpionforcefield.cpp',
               'src/minimizers/lbfgs_interface.cpp',
               'src/minimizers/routines.c',
               'src/minimizers/lbfgs_wrapper/lbfgsb_wrapper.cpp']

    sources.append("bindings/_ptools.pyx")

    ptools = Extension('_ptools',
                       sources=sources,
                       language='c++',
                       library_dirs=[boost_include_dir],
                       include_dirs=['headers',
                                     f2c_include_dir, boost_include_dir],
                       extra_objects=[f2clib])

    cgopt = Extension('cgopt',
                      sources=['PyAttract/cgopt.pyx',
                               'PyAttract/chrg_scorpion.c'],
                      language='c',
                      include_dirs=[f2c_include_dir, 'PyAttract'],
                      extra_objects=[f2clib])

    setup(ext_modules=[ptools, cgopt],
          cmdclass={'build_ext': build_ext},
          packages=['.'],
          name='ptools',
          version='1.2')


if __name__ == '__main__':
    setup_package()
