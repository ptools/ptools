# -*- coding: utf-8 -*-
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import os
import sys
import bzrrev



#### PLEASE override the following two variables
#### to define a custom path to required dependencies:

# user-defined path to libf2c.a
user_path_libf2c = ""

# user-defined path to f2c.h
user_path_f2c_h = ""

# user-defined path to boost library:
user_path_boost = ""



############################




### methods to locate headers and libraries:

def find_file(name,path, useEnvPath=False):
    """finds a file named 'name' in a serie of directories
    given in 'path' """
    
    if useEnvPath:
        additional_path = os.environ['PATH'].split(':')
        path.extend(additional_path)

    for p in path:
        fullfilepath = os.path.join(p,name)
        if os.path.exists(fullfilepath):
            return fullfilepath

    return None #no file found



def find_header(names, paths, useEnvPath=False):
   #find a library in a given set of directories
   
   if useEnvPath:
      additional_path = os.environ['PATH'].split(':')
      paths.extend(additional_path)
   
   for p in paths:
      if os.path.exists(p):
         for n in names:
            if os.path.exists(os.path.join(p,n)):
               return p
                       
   return None



print "Trying to locate the boost libraries, modify setup.py to change de default search path"
boostdir=find_header(["boost/shared_array.hpp"], [user_path_boost,"./ptools_dep/boost_1_55_0", "/usr/include", \
"/opt/local/include"] )

print "Trying to locate libf2c.a static library"
f2clib = find_file("libf2c.a", [user_path_libf2c, "./ptools_dep/f2c", "/usr/lib/", "/usr/local/lib/", "/usr/lib64/"])
f2c_header = find_header(["f2c.h"], [user_path_f2c_h, "./ptools_dep/f2c", "/usr/include", "/usr/local/f2c/"])

not_found_message = """Note: You can add a custom search path by editing this file (setup.py).
You can also install locally the missing dependencies by running: sh ./install-deps.sh"""


if boostdir is None:
    print "Cannot locate boost library", not_found_message
    sys.exit(1)
if f2clib is None:
    print "Cannot locate libf2c", not_found_message
    sys.exit(1)


print "using boost from", boostdir
print "using f2clib from", f2clib



sources =   [
               "cython_wrappers.cpp", 
               "atom.cpp", 
               "attractrigidbody.cpp", 
               "coordsarray.cpp", 
               "mcopff.cpp", 
               "rigidbody.cpp", 
               "surface.cpp",
               "atomselection.cpp", 
               "basetypes.cpp",
               "forcefield.cpp",
               "pairlist.cpp",
               "rmsd.cpp",
               "version.cpp",
               "attractforcefield.cpp",
               "coord3d.cpp",
               "geometry.cpp",
               "pdbio.cpp",
               "superpose.cpp",
               "scorpionforcefield.cpp",
               "imcforcefield.cpp",
               "minimizers/lbfgs_interface.cpp",
               "minimizers/routines.c",
               "minimizers/lbfgs_wrapper/lbfgsb_wrapper.cpp",
               
             ]
             
sources = [os.path.join('src', i) for i in sources ]  #append the 'src' prefix to all files

sources.append("bindings/_ptools.pyx")


setup(ext_modules=[Extension(
                   "_ptools",                 # name of extension
                  sources = sources,
                   language="c++",   # causes Cython to create C++ source
                   #libraries=['f2c'],
                   library_dirs = [boostdir],
                   extra_objects = [f2clib],
                   include_dirs = ['headers', f2c_header, boostdir],
                   )

                     ], 
      cmdclass={'build_ext': build_ext},

      packages = ['.'],
      name="ptools",
      version = "1.2"
      

      
      )
