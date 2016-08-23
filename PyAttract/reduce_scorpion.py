#!/usr/bin/env python
# -*- coding: utf-8 -*-


from optparse import OptionParser
from cgopt import optimize
parser = OptionParser()

# --cgopt option: choice of charge optimization with SCORPION
parser.add_option("--cgopt", dest="optimizecharges", action="store_true", default=False,
                  help="optimize SCORPION coarse grained charges")

parser.add_option("--dgrid", type="float",dest="delgrid", default=1.5,
                  help="grid spacing (A) for charge optimization (default is 1.5), works only with -cgopt option")

# --mcop option: create multicopy reduced protein
parser.add_option("--mcop", dest="regions", default=False, help="positions of the multicopy region(s) separated by ';' (eg. 12-21;45-48)")

(options, args) = parser.parse_args()

import sys
import copy

from ptools import *

try:
    #Python 2.7+
    from collections import OrderedDict
except:
    #for older Python versions with the 'ordereddict' package installed
    from ordereddict import OrderedDict


class IncompleteBead:
      pass

class BeadCreator:

      def __init__(self, reducedname, reducedtypenb, reducedcharge,lstofAtoms):
            self._reducedname=reducedname
            self._reducedtypenb=reducedtypenb
            self._reducedcharge=reducedcharge
            self._lstofAtoms=copy.deepcopy(lstofAtoms)
            self._CoM=Coord3D()  #from ptools
            self.size=0

            atProp=Atomproperty()
            atProp.atomType = reducedname
            atProp.atomCharge = reducedcharge
            atProp.chainId = ''
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
      def create(self):
            "creates a new atom bead"
            if len(self._lstofAtoms)!=0:
                  raise IncompleteBead
            CoM=self._CoM*(1.0/float(self.size))
            at=Atom(self.atProp,CoM)
            return at

###insert load parameter file here


#read parameter file:
import os
parameters = open(os.path.split(sys.argv[0])[0] + "/reduce_data/at2cg_scorpion.dat", 'r').readlines()

beadCorresp = {}

#list of all residues names
#keys are residues names
#values are coarse grain names 
#ie: residNames['ALA'] == {'CA':1 ,'CB':1 }
residNames = {}

#grainMap:
#keys are like this: 'ALA:CA', 
grainMap = {}

initalBeadCharge= {
      'ASP:CB': -1,
      'GLU:CG': -1,
      'ARG:CG': +1, 
      'LYS:CG': 1,
      'GTP:CL': -1,
      'GTP:CM': -1,
      'GTP:CN': -2,
      'GDP:CL': -1,
      'GDP:CM': -2,
      'MG2:CB': +2
}



allAtomRadius = {}
allAtomCharges = {}
beadRadius = {}


for p in parameters:
  lspl = p.split() 
  try:  
    if lspl[0] != '#':
        listofgrains = residNames.get(lspl[0], OrderedDict() )
        listofgrains[lspl[5]] = 1    #add coarse grain name    
        residNames[lspl[0]] = listofgrains

        
        listOfAtomsInAGrain = grainMap.get("%s:%s"%(lspl[0], lspl[5]), [] ) 

        listOfAtomsInAGrain.append(dict(  atomname=lspl[1],  
                                          atomcharge=float(lspl[2]),
                                          atomradius=float(lspl[3]),
                                          weight=float(lspl[4]),
                                          beadname=lspl[5],
                                          beadid=int(lspl[6]),
                                          beadradius=float(lspl[7]),
                                        ) 
                                   )

        allAtomRadius["%s:%s"%(lspl[0], lspl[1]) ] = float(lspl[3]) 
        allAtomCharges["%s:%s"%(lspl[0], lspl[1]) ] = float(lspl[2])
        beadRadius["%s:%s"%(lspl[0], lspl[5])  ] = float(lspl[7]) 

        grainMap["%s:%s"%(lspl[0], lspl[5]) ] = listOfAtomsInAGrain

 
  except IndexError:
    pass



beadCorresp = {}

for residname, cgnames in residNames.items():
    for cgname in cgnames:
        key = "%s:%s"%(residname, cgname)
        lstOfAllAtoms = [i['atomname'] for i in grainMap[key] if i['weight']!=0 ]
        beadid = grainMap[key][0]['beadid']
        beadcharge = initalBeadCharge.get(key, 0.0)
        beadname = grainMap[key][0]['beadname']

        beadDescription = [beadname, lstOfAllAtoms, beadid, beadcharge]
    
        descriptions =  beadCorresp.get(residname, [])
        descriptions.append(beadDescription)
        beadCorresp[residname] = descriptions

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
        write_mcop('_mcop.pdb', core, region_copies)
        models.append(core)
        models.append(region_copies)
else:
    models.append(Rigidbody(args[0]))




def reduce_model(allAtom):  
    newallAtom = []
    for i in xrange(len(allAtom)):
        atom = allAtom.CopyAtom(i)
        if atom.chainId == '': atom.chainId = ' '
        if atom.atomType[0] != 'H' and atom.atomType != 'OXT' and atom.atomType!= 'OT2':
             newallAtom.append(atom)
    allAtom = Rigidbody()
    for at in newallAtom:
        allAtom.AddAtom(at)



    sys.stderr.write("Number of atoms: %d\n" %(len(allAtom) ))

    #extract all 'atoms' objects
    atoms=[]
    for i in xrange(len(allAtom) ):
          atoms.append(allAtom.CopyAtom(i))



    #count residues:
    residuMap={}
    residulist=[]
    for at in atoms:
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

    protein = [] 

    for residKey, atomList in zip(residulist,orderedresid):
          residType=residKey[:3]
          if (residType)=="HIE": residType="HIS" #fix for an amber output file
          residNumber=int(residKey[4:])
          correspList=beadCorresp[residType]

          for correspUnit in correspList:
                atomTypeName=correspUnit[0]
                lstToReduce=correspUnit[1]
                atomTypeNumber=correspUnit[2]
                atomCharge=correspUnit[3]
                beadcreator=BeadCreator(atomTypeName,atomTypeNumber, atomCharge, lstToReduce)
                for atom in atomList:
                      beadcreator.submit(atom)
                try:
                      bead = beadcreator.create()
                except IncompleteBead:
                      sys.stderr.write("The bead %i of residue %s is incomplete. Please check your pdb!\n"\
                          %(totAtoms+1,residKey) )
                      sys.exit(1)
                totAtoms+=1
                #now we must modify the bead: change the residue type and set the "extra" field correctly
                bead.residType = residType 
                extra = ('%5i'+'%8.3f'+'%2i'*2) %(atomTypeNumber,atomCharge,0, 0)
                bead.extra = extra
                bead.atomId = totAtoms
                bead.residId = residNumber
                protein.append(bead)


    charge = []
    radius = []
    cx = []
    cy = []
    cz = []

    cgch = []
    cgr = []
    cgx = []
    cgy = []
    cgz = []



    for i in range(len(allAtom)):
       atom = allAtom.CopyAtom(i)
       residu_type= atom.residType
       atomtype = atom.atomType
       if residu_type =="ILE" and atomtype == "CD":
           atomtype = "CD1"
       key = "%s:%s"%(residu_type, atomtype) 
       radius.append(  allAtomRadius[key] ) 
       charge.append ( allAtomCharges[key]   )
       cx.append( atom.coords.x)
       cy.append( atom.coords.y)
       cz.append( atom.coords.z)

       
      
    for i, atom in enumerate(protein):
       
       cgch.append( atom.atomCharge ) 
       
       residu_type= atom.residType
       atomtype = atom.atomType
       if residu_type =="ILE" and atomtype == "CD":
           atomtype = "CD1"
       key = "%s:%s"%(residu_type, atomtype) 

       cgr.append( beadRadius[key] )  
       cgx.append( atom.coords.x)
       cgy.append( atom.coords.y)
       cgz.append( atom.coords.z)






    first = False
    last = False
    for i, at in enumerate(protein):
       if at.atomType == 'CA':
          cgch[i] += 1
          break
    protein.reverse()
    cgch.reverse()

    for i,at in enumerate(protein):
       if at.atomType == 'CA':
          cgch[i] -= 1
          break
          
    protein.reverse()
    cgch.reverse()  


    if options.optimizecharges:
        optimized = optimize(len(allAtom), charge, radius, cx, cy, cz, len(protein), cgch, cgr, cgx, cgy, cgz, options.delgrid )
    else:
        optimized = cgch


    for i, bead in enumerate(protein):
        #ugly hack to set correct charges values due to a bug atom to pdb conversion
        extra = bead.extra
        atomTypeNumber = int(extra.split()[0])
        bead.extra = ('%5i'+'%8.3f'+'%2i'*2) %(atomTypeNumber,optimized[i],0, 0)
        
        print bead.ToPdbString()



# iterate through all models (important for multicopy option) and write reduced pdb
print "HEADER    SCORPION REDUCED PDB FILE"
# iterate through 'models' variable which contains max 2 elements.
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
            
        
        
