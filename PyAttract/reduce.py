#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import copy
import os
from optparse import OptionParser

usage = "%prog --ff force_field [--dna] [--cgopt] [--dgrid 1.5] [--allow_missing] atomic_file.pdb > reduced_file.red"
version = "%prog 1.0"
parser = OptionParser(usage)

# --ff option: choice of forcefield
parser.add_option("--ff",type="str",dest="ffversion",help="choice of CG force field: attract1, attract2, scorpion, or imc")

# --cgopt option: choice of charge optimization with SCORPION
parser.add_option("--cgopt", dest="optimizecharges", action="store_true", default=False,
                  help="optimize SCORPION coarse grained charges, works only with FF Scorpion")

parser.add_option("--dgrid", type="float",dest="delgrid", default=1.5,
                  help="grid spacing (A) for charge optimization (default is 1.5), works only with FF Scorpion and -cgopt option")

# --dna option: reduce dna
parser.add_option("--dna", action="store_true", dest="molDna", default=False, help="reduce dna, works only with FF Attract1")

parser.add_option("--allow_missing", action="store_true", dest="warning",default=False, help="don't stop program if atoms are missing, only display a warning on stderr")


(options, args) = parser.parse_args()

#==========================================================
# check options
#==========================================================

#get the location of this script:
thisscript=os.path.realpath(__file__)
thispath = os.path.split(thisscript)[0]
print "thispath is: ", thispath

# define other parameter files
cmd_options=[]
if options.ffversion not in ["attract1","attract2","scorpion","imc"]:
    sys.stderr.write ("Error: please choose one of the following CG force field: attract1, attract2, scorpion, or imc\n")
    sys.exit(1)
if options.warning:
    if options.ffversion!="scorpion":
       cmd_options.append("--allow_missing")
    else:
       sys.stderr.write("Warning: allow_missing option not supported by SCORPION force field\n")
if options.ffversion=="scorpion":
    if options.optimizecharges:
       cmd_options.append("--cgopt")
    if options.delgrid:
       cmd_options.append("--dgrid %f"%options.delgrid)
    cmd_options.append(args[0])
    prgname = os.path.join(thispath, "reduce_scorpion.py")

    os.system(prgname + " " + " ".join(cmd_options))

if options.ffversion=="attract1":
    if options.molDna:
       cmd_options.append("--dna")
    else:
       cmd_options.append("--prot")
    cmd_options.append(args[0])
    prgname = os.path.join(thispath, "reduce_attract1.py")

    os.system(prgname +  " %s"%" ".join(cmd_options))

if options.ffversion=="attract2":
    cmd_options.append(args[0])
    prgname = os.path.join(thispath, "reduce_attract2.py")

    os.system(prgname +  " %s"%" ".join(cmd_options))

if options.ffversion=="imc":
    cmd_options.append(args[0])
    prgname = os.path.join(thispath, "reduce_imc.py")

    os.system(prgname +  " %s"%" ".join(cmd_options))
