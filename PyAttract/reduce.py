#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import copy
import os
from optparse import OptionParser

# ---- Functions ---- #
def is_digit_array(s):
    for digit in s:
        if not digit.isdigit():
            return False
    return True
# -- End of Functions -- #

usage = "%prog --ff force_field [--dna] [--cgopt] [--dgrid 1.5] [--allow_missing] atomic_file.pdb > reduced_file.red"
version = "%prog 1.0"
parser = OptionParser(usage)

# --ff option: choice of forcefield
parser.add_option("--ff",type="str",dest="ffversion",help="choice of CG force field: attract1, attract2 or scorpion")

# --cgopt option: choice of charge optimization with SCORPION
parser.add_option("--cgopt", dest="optimizecharges", action="store_true", default=False,
                  help="optimize SCORPION coarse grained charges, works only with FF Scorpion")

parser.add_option("--dgrid", type="float",dest="delgrid", default=1.5,
                  help="grid spacing (A) for charge optimization (default is 1.5), works only with FF Scorpion and -cgopt option")

# --dna option: reduce dna
parser.add_option("--dna", action="store_true", dest="molDna", default=False, help="reduce dna, works only with FF Attract1")

parser.add_option("--allow_missing", action="store_true", dest="warning",default=False, help="don't stop program if atoms are missing, only display a warning on stderr")

# --mcop option: create multicopy reduced protein
parser.add_option("--mcop", dest="regions", default=False, help="positions of the multicopy region(s) separated by ':' (eg. --mcop 12-21:45-48)\n If entire protein, enter ':' (eg. --mcop :")


(options, args) = parser.parse_args()

#==========================================================
# check options
#==========================================================

#get the location of this script:
thisscript=os.path.realpath(__file__)
thispath = os.path.split(thisscript)[0]
print "thispath is: ", thispath

# define other parameter files

# Check force field option (--ff)
cmd_options=[]
if options.ffversion not in ["attract1","attract2","scorpion"]:
    sys.stderr.write ("Error: please choose one of the following CG force field: attract1, attract2 or scorpion\n")
    sys.exit(1)
if options.warning:
    if options.ffversion!="scorpion":
       cmd_options.append("--allow_missing")
    else:
       sys.stderr.write("Warning: allow_missing option not supported by SCORPION force field\n")


# Check multicopy option (--mcop)
if options.regions:
    if options.regions.endswith('.pdb'):
        cmd_options.append('--mcop .pdb')
        args.insert(0,options.regions)
    elif options.regions == ':':
        cmd_options.append("--mcop :")
    else:
        pos = options.regions.replace(':','-').split('-')
        reg = options.regions.split(':')
        num_positions = len(pos)
        num_regions = len(reg)

        if num_positions != 2*num_regions or not is_digit_array(pos):
            sys.stderr.write ("Error: the --mcop input format is incorrect. \n Separate start and end of a region with '-'. \n Separate different regions with ':' \n (Example for 2 regions : --mcop 12-21:45-48)\n Or enter the files containing core, and region copies.\n (Example for 2 regions : --mcop core.pdb region1.pdb region2.pdb)\n")
            sys.exit(1)

        starts = []
        ends = []
        for i in xrange(0,len(reg)):
            starts.append(min(map(int, reg[i].split('-'))))
            ends.append(max(map(int, reg[i].split('-'))))
        starts, ends = zip(*sorted(zip(starts,ends)))

        for i in xrange(0,len(reg)-1):
            if ends[i] >= starts[i+1]:
                sys.stderr.write("Error: the --mcop input has overlapping regions.\n")
                sys.exit(1)

        cmd_options.append("--mcop %s"%options.regions)

#Launch the selected force field version. 
if options.ffversion=="scorpion":
    if options.optimizecharges:
       cmd_options.append("--cgopt")
    if options.delgrid:
       cmd_options.append("--dgrid %f"%options.delgrid)
    cmd_options.extend(args)
    prgname = os.path.join(thispath, "reduce_scorpion.py")

    os.system(prgname + " " + " ".join(cmd_options))

if options.ffversion=="attract1":
    if options.molDna:
       cmd_options.append("--dna")
    else:
       cmd_options.append("--prot")
    cmd_options.extend(args)
    prgname = os.path.join(thispath, "reduce_attract1.py")

    os.system(prgname +  " %s"%" ".join(cmd_options))

if options.ffversion=="attract2":
    cmd_options.extend(args)
    prgname = os.path.join(thispath, "reduce_attract2.py")

    os.system(prgname +  " %s"%" ".join(cmd_options))
