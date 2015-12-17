#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import copy
import os
from optparse import OptionParser


from ptools import *


usage = "%prog --prot --dna atomic_file.pdb [--red file] [--ff file] [--conv file] [--allow_missing] > reduced_file.red"
version = "%prog 1.0"
parser = OptionParser(usage)

# --prot option: reduce protein
parser.add_option("--prot", action="store_true", dest="molProt", default=False, help="reduce protein")

# --dna option: reduce dna
parser.add_option("--dna", action="store_true", dest="molDna", default=False, help="reduce dna")

# --red option: correspondance between atoms and beads (coarse grain)
parser.add_option("--red", action="store", type="string", dest="redName", help="correspondance file between atoms and beads")

# --ff option: coarse grain force field parameters
#parser.add_option("--ff", action ="store", type="string", dest="ffName", help="force field parameters") # imc

## --conv option: eventually, a correspondance between different names for the same atom
#parser.add_option("--conv", action='store', type='string', dest='convName', help='type conversion file') # imc

parser.add_option("--allow_missing", action="store_true", dest="warning",default=False, help="don't stop program if atoms are missing, only display a warning on stderr")



(options, args) = parser.parse_args()

#==========================================================
# check options
#==========================================================

# define reduce data files subdirectory
data_dir="reduce_data/"

if ( (options.molProt or options.molDna) and ( len(args) > 0) ):
        atomicName = args[0]
        completePath=sys.argv[0]
        scriptdir,scriptname = os.path.split(completePath)
        if options.molProt:
                redName = os.path.join(scriptdir,data_dir+"at2cg_imc.dat")
        if options.molDna:
                redName = os.path.join(scriptdir,data_dir+"at2cg.dna.dat")
        ffName = os.path.join(scriptdir,data_dir+"ff_param.dat") # imc
        #convName = os.path.join(scriptdir,data_dir+"type_conversion.dat") # imc  
else:
        parser.error("please specify molecule type (--prot or --dna) and atomic file")

# define other parameter files
if options.redName:
        redName  = options.redName
#if options.ffName: # imc
#        ffName   = options.ffName # imc
#if options.convName: # imc
#        convName = options.convName # imc

#==========================================================
# check files
#==========================================================
# check if a required file is found
def checkFile(name):
        flag = os.path.exists(name)
        if  not flag :
                sys.stderr.write("ERROR: cannot find %s \n" %(name))
                exit(2)	
checkFile(atomicName)
checkFile(redName)
#checkFile(ffName) # imc
#checkFile(convName) # imc



#==========================================================
# classe Atom
#==========================================================
class AtomInBead:
        'Class definition for an atom in a bead'
        def __init__(self, name, wgt):
                self.name = name    # atom name
                self.x = 0.0        # x coordinate
                self.y = 0.0        # y coordiante
                self.z = 0.0        # z coordinate
                self.weight = wgt   # atom weight (within a bead)
                self.found = 0      # atom found or not (not found by default)

        def Show(self):
                return "%s %8.3f %8.3f %8.3f %3.1f %d" %(self.name, self.x, self.y, self.z, self.weight, self.found)


#==========================================================
# class Bead
#==========================================================
class Bead:
        'Class definition for a bead'
        def __init__(self, name, id):
                self.name = name        # bead name
                self.id = id            # bead id
                self.size = 0           # bead size (number of atoms inside)
                self.listOfAtomNames = []  # list of all atom names
                self.listOfAtoms = []  # list of all atoms (AtomInBead)

        def Show(self):
                return "%s %d %d" %(self.name, self.id, self.size)


#==========================================================
# class Res
#==========================================================
class CoarseRes:
        'Class definition for coarse grain (reduced) protein residue (or DNA base)'
        def __init__(self):
                self.listOfBeadId = []   # list of bead id in res
                self.listOfBeads = []       # list of beads in res

        def Add(self, residue):
                'Add in a residue atoms from bead'
                for at in residue:
                        at_name = at[1]
                        at_wgt = float(at[2])
                        bd_id = int(at[3])
                        bd_name = at[4]
                        if at_name != 'EMPTY': # EMPTY is a special tag to deal with glycine
                                # in bead not in residue than create it
                                if bd_id not in self.listOfBeadId:
                                        self.listOfBeadId.append(bd_id)
                                        self.listOfBeads.append(Bead(bd_name, bd_id))
                                # add atom in bead in residue
                                bead_position = self.listOfBeadId.index(bd_id)
                                bead = self.listOfBeads[bead_position]
                                bead.listOfAtomNames.append(at_name)
                                atInBd = AtomInBead(at_name, at_wgt)
                                bead.listOfAtoms.append(atInBd)
                                bead.size += 1
                                # update bead in residue
                                self.listOfBeads[bead_position] = bead
                # return the number of bead per residue
                return 	len(self.listOfBeadId)	


        def FillAtom(self, at_name, x, y, z):
                'Fill an atom from bead with coordinates'
                # quickly check atom in atom list
                # 1: browse beads
                for bead in self.listOfBeads:
                        # 2: browse atoms in bead
                        if at_name in bead.listOfAtomNames:
                                # then find exactly where this atom is present
                                for atom in bead.listOfAtoms:
                                        if at_name == atom.name:
                                                atom.x = x
                                                atom.y = y
                                                atom.z = z
                                                atom.found = 1
                                                #print "fill", atom.name, atom.x, atom.y, atom.z


        def Reduce(self, infoResName, infoResId):
                'Reduce a bead with atoms present in bead'
                output = []
                # reduce all beads in a residue
                # for each bead in the residue
                for bead in self.listOfBeads:
                        #print "reducing", bead.name, bead.id, "with", bead.listOfAtomNames
                        reduce_size = 0
                        reduce_x = 0.0
                        reduce_y = 0.0
                        reduce_z = 0.0
                        sum_wgt = 0.0
                        # for each atom of a bead
                        for atom in bead.listOfAtoms:
                                if atom.found == 1:
                                        reduce_size += 1
                                        reduce_x += atom.x * atom.weight
                                        reduce_y += atom.y * atom.weight
                                        reduce_z += atom.z * atom.weight
                                        sum_wgt += atom.weight
                                else:
                                        message="ERROR: missing atom %s in bead %s %2d for residue %s %d. Please fix your PDB!\n" \
                                        %(atom.name, bead.name, bead.id, infoResName, infoResId)
                                        if options.warning:
                                            sys.stderr.write(message)
                                            sys.stderr.write("Continue execution as required ...\n")
                                        else:
                                            raise Exception(message)
                        if reduce_size == bead.size:
                                coord = Coord3D(reduce_x/sum_wgt, reduce_y/sum_wgt, reduce_z/sum_wgt)
                                output.append([coord, bead.name, bead.id])				
                return output

        def Show(self):
                for bead in self.listOfBeads:
                        print bead.name, bead.id, bead.size, bead.atomIdList

#==========================================================
# definition of the ATTRACT coarse grain model
# needs beads definition from redName file 
# and bead charges from ff_name file
#==========================================================
f = open(redName, 'r')
lines = f.readlines()
f.close()
#
# ARG       CG        1.0       3         CG
#
resBeadAtomModel = {}

# quickly parse file to get full list of protein residues / DNA bases
# and to create an empty dictionnary with residues or bases as keys
sys.stderr.write("%s: found the definition of residues " %(redName))
for line in lines:
        if (len(line) > 1) and (line[0] != '#'):
                item = line.split()
                if len(item) >= 5: # at least 5 fields are expected
                        res = item[0]
                        if (res != '*') and (res not in resBeadAtomModel.keys()):
                                resBeadAtomModel[res] = CoarseRes()
                                sys.stderr.write('%s ' %(res))
sys.stderr.write('\n')

# fill for each residue or base, the beads and atoms
sys.stderr.write("%s: created the partition for residues " %(redName))
beadsInResidue = []
for line in lines:
        if (len(line) > 1) and (line[0] != '#'):
                col = line.split()
                if len(col) >= 5: # at least 5 fields are expected
                        beadsInResidue.append(col[0:5])
        if (line[0:4] == '#===') and (len(beadsInResidue) != 0):
                if beadsInResidue[0][0] == "*":
                        # "*" means a common bead for all residues or bases
                        for res in resBeadAtomModel.keys():
                                resBeadAtomModel[res].Add(beadsInResidue)
                else:
                        # if not a common bead, add it just to the right residue (or base)
                        res = beadsInResidue[0][0]
                        bead_nb = resBeadAtomModel[res].Add(beadsInResidue)
                        sys.stderr.write('%s(%d beads) ' %(res, bead_nb))
                beadsInResidue = []
sys.stderr.write('\n')


#==========================================================
# read force field parameter file to get bead charge
#==========================================================
#f = open(ffName, 'r')
#lines = f.readlines()
#f.close()
#sys.stderr.write('%s: reading force field parameters for bead ' %(ffName))
beadChargeDic = {}
#for line in lines:
#        if (len(line) > 1) and (line[0] != '#'):
#                item = line.split()
#                if len(item) < 5: # at least 5 fields are expected
#                        break
#                beadId = int(item[0])
#                beadCharge = float(item[3])
#                beadChargeDic[beadId] = beadCharge
#                sys.stderr.write('%d ' %(beadId))
#sys.stderr.write('\n')


#==========================================================
# read file for type conversions (residues or atoms)
#==========================================================
resConv = {}
atomConv = {}
#f = open(convName, 'r')
#lines = f.readlines()
#f.close()
#for line in lines:
#        if (len(line) > 1) and (line[0] != '#'):
#                item = line.split()
#                if len(item) == 2: # 2 fields => residue type conversion
#                        resOld = item[0]
#                        resNew = item[1]
#                        if resOld not in resConv.keys():
#                                resConv[resOld] = resNew
#                if len(item) == 3: # 3 fields => atom type conversion
#                        res = item[0]
#                        atomOld = item[1]
#                        atomNew = item[2]
#                        atomTagOld = res + '-' + atomOld
#                        atomTagNew = res + '-' + atomNew
#                        if atomTagOld not in atomConv.keys():
#                                atomConv[atomTagOld] = atomTagNew


#==========================================================
# load atomic pdb file into Rigidbody object
#==========================================================
allAtom=Rigidbody(atomicName)
sys.stderr.write("Load atomic file %s with %d atoms \n" %(atomicName, len(allAtom)))

#extract all 'atoms' objects
atomList=[]
for i in xrange(len(allAtom)):
        atom = allAtom.CopyAtom(i)
        # look for residue or base type conversion
        resName = atom.residType
        if resName in resConv.keys():
                atom.residType =  resConv[resName] 
        # look for atom type conversion
        atomTag = atom.residType + '-' + atom.atomType
        if atomTag in atomConv.keys():
                atomName = atomConv[atomTag].split('-')[1] 
                atom.atomType =  atomName 
        atomList.append(atom)

#count residues
residueTagList=[]
coarseResList=[]
for atom in atomList:
        resName = atom.residType
        # create a unique identifier for every residue
        # resTag is for instance "LEU-296-A"
        resTag = resName + '-'+ str(atom.residId) + '-' + atom.chainId
        if resTag not in residueTagList:
                if resBeadAtomModel.has_key(resName):
                        residueTagList.append(resTag)
                        # add a pattern residue to the list of coarse residues for the protein
                        # beware of the hugly list copy: use copy.deepcopy() !
                        coarseResList.append(copy.deepcopy(resBeadAtomModel[resName]))
                else:
                        sys.stderr.write("WARNING: residue %s is unknown the residues <-> beads <-> atoms list !!\n" %(resName))
                        sys.stderr.write("       : residue %s will not be reduced into coarse grain\n" %(resName))
sys.stderr.write("Number of residues: %i\n" %(len(residueTagList)))

#==========================================================
# iterate through all atoms and residues to fill beads
#==========================================================
sys.stderr.write("Reading all atoms and filling beads:\n")
for atom in atomList:
        #resTag is like "LEU-296-A"
        resTag = atom.residType + '-' + str(atom.residId) + '-' + atom.chainId
        if resTag in residueTagList:
                id = residueTagList.index(resTag)
                coarseResList[id].FillAtom(atom.atomType, atom.coords.x, atom.coords.y, atom.coords.z)
#==========================================================
# reduce beads
#==========================================================
coarsegrainPdb = ""   # complete coarse grain (reduced) pdb file
atomCnt = 0           # atom counter
sys.stderr.write("Coarse graining:\n")
for i in range(len(residueTagList)):
        tag = residueTagList[i].split('-')
        resName = tag[0]
        resId = int(tag[1])
        coarseRes = coarseResList[i].Reduce(resName, resId)
        for bead in coarseRes:
                coord = bead[0]
                atomName = bead[1]
                atomTypeId = bead[2]
                if atomTypeId in beadChargeDic:
                        atomCharge = beadChargeDic[atomTypeId]
                else:
                        #sys.stderr.write("WARNING: cannot find charge of bead %s %2d in %s \n" %(atomName, atomTypeId, ffName)) # imc
                        #sys.stderr.write("       : set default charge to 0.0\n") # imc
                        atomCharge = 0.0
                prop = Atomproperty()
                prop.atomType = atomName
                atomCnt += 1
                prop.atomId = atomCnt
                prop.residId = resId
                prop.residType = resName
                prop.chainId = ' '
                extra = ('%5i%8.3f%2i%2i') %(atomTypeId,atomCharge,0,0)
                prop.extra = extra
                newAtom = Atom(prop, coord)
                coarsegrainPdb += newAtom.ToPdbString() + "\n"

#==========================================================
# output coarse grain (reduced) pdb file
#==========================================================
sys.stdout.write("HEADER    IMC REDUCED PDB FILE\n")
sys.stdout.write(coarsegrainPdb)
sys.stderr.write("Coarse grain (reduced) output")
sys.stderr.write(": %d beads \n" %(atomCnt))
