#!/usr/bin/env python
# -*- coding: utf-8 -*-

#modified reduce_ff2 for mbest1u 
#Calpha are dummy atom type 32. 
#Gly Calpha have dummy atom type 1

import sys
import copy

from ptools import *



from optparse import OptionParser
parser = OptionParser()
parser.add_option("--allow_missing", dest="ignoremissing", action="store_true", default=False,
                  help="ignore missing heavy atoms (which will result in missing beads)" )

# --mcop option: create multicopy reduced protein
parser.add_option("--mcop", dest="regions", default=False, help="positions of the multicopy region(s) separated by ';' (eg. 12-21;45-48)")

(options, args) = parser.parse_args()




class IncompleteBead:
      pass

class BeadCreator:

      def __init__(self, reducedname, reducedtypenb, reducedcharge,lstofAtoms, chainId=''):
            self._reducedname=reducedname
            self._reducedtypenb=reducedtypenb
            self._reducedcharge=reducedcharge
            self._lstofAtoms=copy.deepcopy(lstofAtoms)
            self._CoM=Coord3D()  #from ptools
            self.size=0

            atProp=Atomproperty()
            atProp.atomType = reducedname
            atProp.atomCharge = reducedcharge
            atProp.chainId = chainId
            self.atProp=atProp

      def submit(self, atom):
            "try to add an atom to the bead"
            atomtype=atom.atomType
            #trick to handle 'OTn' instead of 'O' for last pdb atom:
            if atomtype[:2]=='OT': atomtype='O'
            
            if atom.residType =="ILE" and atomtype == "CD":
                atomtype = "CD1"
                atom.atomType = "CD1"
            
            if atomtype in self._lstofAtoms:
                  self._CoM+=atom.coords
                  self._lstofAtoms.remove(atomtype)
                  self.size+=1
                  #print self.size
      def create(self):
            "creates a new atom bead"
            if len(self._lstofAtoms)!=0:
                  raise IncompleteBead
            CoM=self._CoM*(1.0/float(self.size))
            at=Atom(self.atProp,CoM)
            return at



defaultBB=[ [ 'N',    ['N']  , 30 , 0.0 ],
            ['CA',    ['CA'] ,  32 , 0.0 ],
            ['C' ,    ['C']  ,  32, 0.0  ],
            ['O',     ['O']  ,  31, 0.0  ]
           ]

beadCorresp={}
beadCorresp["ARG"] = defaultBB + [[  'CG',        ['CG']        ,   3 ,  0.0  ],
                                 [  'NEC',      ['NE','CZ']     ,   4 ,  1.0   ]]

beadCorresp["GLU"] = defaultBB + [
                                  [   'CB',       ['CG']        ,   10  , 0.0  ] ,
                                  [ 'CO1',  ['CD', 'OE1', 'OE2'],   11  , -1.0 ]
                                  ]

beadCorresp["GLN"] = defaultBB + [
                                   [ 'CB',        ['CG']        ,   8,    0.0 ],
                                   [ 'CN1', ['CD', 'OE1', 'NE2'],   9,    0.0 ]
                                  ]


beadCorresp["LYS"] = defaultBB + [
                                   ['CB', ['CG'], 16, 0.0],
                                   ['CE', ['CE'], 17, 1.0]
                                 ]

beadCorresp["TRP"] = defaultBB + [
                                   ['CG', ['CG'], 25, 0.0],
                                   ['CSE', ['CD2','CE2','CE3','CH2','CZ3','CZ2'], 26, 0.0]
                                 ]

beadCorresp["MET"] =  defaultBB + [
                                    ['CSE', ['CB','CG'], 18, 0.0],
                                    ['CSE', ['SD','CE'], 19, 0.0]
                                  ]

beadCorresp["PHE"] = defaultBB + [
                                   ['CSE', ['CB','CG'], 20, 0.0],
                                   ['CSE', ['CD1','CD2','CE1','CE2','CZ'], 21, 0.0 ]
                                 ]

beadCorresp["TYR"] = defaultBB + [  ['CSE', ['CB','CG'], 27, 0.0],
                                    ['CSE', ['CD1','CD2','CE1','CE2','CZ','OH'],28,0.0]  ]

beadCorresp["HIS"] = defaultBB + [  ['CSE', ['CB','CG'], 12, 0.0],
                                    ['CSE', ['ND1','CD2','NE2','CE1'], 13, 0.0]  ]

beadCorresp["GLY"] = copy.deepcopy(defaultBB)
assert(beadCorresp["GLY"][1][2] == 32)
beadCorresp["GLY"][1][2] = 1




beadCorresp['ASN'] = defaultBB + [  ['CSE', ['CB', 'CG', 'OD1', 'ND2'] , 5, 0.0]  ]
beadCorresp["ALA"] = defaultBB + [  ['CSE',      ['CB'],                 2, 0.0]]
beadCorresp['ASP'] = defaultBB + [  ['CSE', ['CB', 'CG', 'OD1', 'OD2'] , 6 ,-1.0] ]
beadCorresp['CYS'] = defaultBB + [  ['CSE',    ['CB', 'SG']            , 7 , 0.]  ]
beadCorresp['ILE'] = defaultBB + [  ['CSE', ['CB', 'CG1', 'CG2', 'CD1'], 14, 0.]  ]
beadCorresp['LEU'] = defaultBB + [  ['CSE', ['CB', 'CG', 'CD1', 'CD2'] , 15, 0.]  ]
beadCorresp['PRO'] = defaultBB + [  ['CSE',     ['CB', 'CG', 'CD']     , 22, 0.]  ]
beadCorresp['SER'] = defaultBB + [  ['CSE',     ['CB', 'OG']           , 23, 0.]  ]
beadCorresp['THR'] = defaultBB + [  ['CSE',     ['CB', 'OG1', 'CG2']   , 24, 0.]  ]
beadCorresp['VAL'] = defaultBB + [  ['CSE',     ['CB', 'CG1', 'CG2']   , 29, 0.]  ]

#fix for charmm: add type HSE HSP and HSP
beadCorresp['HSE'] = beadCorresp['HIS']
beadCorresp['HSD'] = beadCorresp['HIS']
beadCorresp['HSP'] = beadCorresp['HIS']



models = []
if options.regions:
# Preparing pdb multicopy models
    if options.regions.endswith('.pdb'):
        core = Rigidbody(args[0])
        region_copies = []
        for i in xrange(1, len(args)):
            region_copies.append(Mcop(args[i]))
        models.append(core)
        models.append(region_copies)
        write_mcop('new_mcop.pdb', core, region_copies)
    else:
        core, region_copies = multicopy(args[0], options.regions)
        write_mcop('new_mcop.pdb', core, region_copies)
        models.append(core)
        models.append(region_copies)
else:
    models.append(Rigidbody(sys.argv[1]))
    sys.stderr.write("Number of atoms:%i \n" %len(models[0]) )

def reduce_model(allAtom):
    #extract all 'atoms' objects
    atoms=[]
    for i in xrange(len(allAtom)):
          atoms.append(allAtom.CopyAtom(i))



    #count residues:
    residuMap={}
    residulist=[]

    # chain id for the reduced file:
    outChainId = atoms[0].chainId

    for at in atoms:
          #fix for incorrect pdb: append a chainId when it's missing
          if at.chainId=='':
              at.chainId = 'A'
          residueIdentifier = at.residType + str(at.chainId)  + str(at.residId)
          #residueIdentifier is like "LEUA296"
          residuMap.setdefault(residueIdentifier, []).append(at)
          if residueIdentifier not in residulist:
                residulist.append(residueIdentifier)

    sys.stderr.write("Number of residues: %i\n" %(len(residuMap)))
    sys.stderr.write("Start atom of each residue:\n")
    orderedresid=[residuMap[i] for i in residulist ]
    startatoms=[lat[0].atomId for lat in orderedresid ]
    out = ""
    for statom in startatoms:
          out+=(str(statom))+" "
    sys.stderr.write(out+"\n")


    #iterates through all the residues and create reduced beads:

    totAtoms=0

    index = 0
    for residKey, atomList in zip(residulist,orderedresid):
          residType=residKey[:3]
          if (residType)=="HIE": residType="HIS" #fix for an amber output file
          #print "key:", residKey
          residNumber=int(residKey[4:])
          #print "reducing residue %s of type %s" % (residKey, residType)
          correspList=beadCorresp[residType]
          #print correspList
          for correspUnit in correspList:
                atomTypeName=correspUnit[0]
                lstToReduce=correspUnit[1]
                atomTypeNumber=correspUnit[2]
                atomCharge=correspUnit[3]
                beadcreator=BeadCreator(atomTypeName,atomTypeNumber, atomCharge, lstToReduce, outChainId )
                for atom in atomList:
                      beadcreator.submit(atom)
                try:
                      bead = beadcreator.create()
                except IncompleteBead:
                      sys.stderr.write("The bead %i of residue %s is incomplete. Please check your pdb!\n"\
                          %(totAtoms+1,residKey) )
                      if not options.ignoremissing: raise 
                totAtoms+=1
                #now we must modify the bead: change the residue type and set the "extra" field correctly
                bead.residType = residType 
                extra = ('%5i'+'%8.3f'+'%2i'*2) %(atomTypeNumber,atomCharge,0, 0)
                bead.extra = extra
                bead.atomId = totAtoms
                bead.residId = residNumber
                print bead.ToPdbString()
                
# iterate through all models (important for multicopy option) and write reduced pdb
print "HEADER    ATTRACT2 REDUCED PDB FILE"
# iterate through 'models' varialbe which contains max 2 elements.
# The fist element of the models list is a Rigidbody (the core in the mcop case, and the whole protein in the normal case)
# In the mcop case, there is a second element, which is the region_copies list (a list of Mcop objects).
for i in xrange(0,len(models)):
    if i == 0 and options.regions != ':':
        sys.stderr.write("\nCore:\n")
        print 'MODEL       0'
        reduce_model(models[i])
        print 'ENDMDL      0'
    if i >= 1:
        # iterate through number of regions
        for j in xrange(0, len(models[1])):
            # iterate through number of regions copies
            for k in xrange(0,len(models[1][j])):
                sys.stderr.write("\nRegion %i Copy %i:\n" %(j+1, k+1))
                print 'MODEL       ' + str(j+1) + '  ' + str(k+1)
                reduce_model(models[1][j][k])
                print 'ENDMDL    ' + str(j+1) + '  ' + str(k+1)
        sys.stderr.write("\n\nSummary:\n")
        for j in xrange(0, len(models[1])):
            sys.stderr.write("Region %i:    %i copies\n" %(j+1,len(models[1][j])))

