#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ptools import *
import sys
import os
#import time
import datetime
import math
import string
import bz2  #for compression of Ligand and receptor data
import base64 #compressed ligand and receptor as base64 strings



def surreal(i):
    return i

def rmsdca(l1,l2):
    return Rmsd(l1.CA().CreateRigid(), l2.CA().CreateRigid())


def compress_file(filename):
    fobject = open(filename,"r")
    all = fobject.read()
    compressed = bz2.compress(all)
    encoded = base64.b64encode(compressed)
    return "compressed %s : \"%s\""%(filename,encoded)


def PrintVect(vect):
    for i in range(len(vect)):
        print vect[i], " | ",
    print ''



class Rotation:
    
    class _Rot:
        ssi=0.0
        phi=0.0
        rot=0.0
        def __init__(self, ssii, phii, roti):
            self.ssi=ssii
            self.phi=phii
            self.rot=roti
    
    def read_rotdat(self):
        self.zwopi=2.0*3.14159265
        self.NbRotByTrans=0
        self.theta=[]
        self.nphi=[]
        # read theta,phi,rot data
        rotdat=open('rotation.dat','r')
        line=rotdat.readline().split()
        self.ntheta=int(line[0])
        self.nrot=int(line[1])
        print "ntheta, nrot: %i %i" %(self.ntheta,self.nrot)
        for i in range(self.ntheta):
            line=rotdat.readline().split()
            self.theta.append(float(line[0]))
            self.nphi.append(int(line[1]))
            self.NbRotByTrans=self.NbRotByTrans+self.nphi[i]*self.nrot
            self.theta[i]=self.zwopi*self.theta[i]/360.0
            print self.theta[i],self.nphi[i]
        rotdat.close()
        self._rot=[]
        
        print "%i rotations by translation"%self.NbRotByTrans

        for kkk in range(self.ntheta):
            ssii=self.theta[kkk]
            phii=self.zwopi/self.nphi[kkk]
            for jjj in range(self.nphi[kkk]):
                phiii=(jjj+1)*phii
                for iii in range(self.nrot):
                    roti=(iii+1)*self.zwopi/self.nrot
                    self._rot.append((phiii, ssii, roti))


    def __init__(self):
        self.read_rotdat()
    
    def __iter__(self):
        return self._rot.__iter__()



class Translation:
    def __init__(self):
        self.translation_dat=Rigidbody("translation.dat")
        print "Reading %i translations from translation.dat"%len(self.translation_dat)

    def __iter__(self):
        self.i=0
        return self
    def next(self):
        if (self.i == len(self.translation_dat)): raise StopIteration
        coord=self.translation_dat.getCoords(self.i)
        self.i+=1
        return [self.i,coord]
        


def skipComments(attach):
    result=[]
    for line in attach.readlines():
        if line[0] not in ["#","!"]:
            result.append(line)
    return result



def readParams(filename):
    fich=open(filename,"r")
    clean=skipComments(fich)
    
    #nb of minimizations to perform is first integer of cleaned lines
    nbminim=int(clean.pop(0).split()[0])
    print "%i series of minimizations"%nbminim
    
    
    lignames=[]
    while(True):
        line=clean.pop(0).split()
        if line[0]=="Lig":
            lignames.append(line[2])
        else:
            break

    rstk = float(line[3])
    #print "rstk = ",rstk
    #ignored=line #one more line is read in the loop before.
                 # (the xrst line which is ignored now)
    
    minimlist=[]
    for i in range(nbminim):
        line=clean.pop(0).split()
        minim={}
        minim['maxiter'] = int(line[0])
        minim['squarecutoff'] = float(line[-1])
        if int(line[-2])==1:
            minim['rstk'] = rstk
        else:
            minim['rstk'] = 0.0
        minimlist.append(minim)  #return a list of type (iter,squarecutoff,has_restraint)
                                 #other parameters are ignored
    return (nbminim, lignames, minimlist,rstk)




def rigidXstd_vector(rigid, mat_std):
    
    #create a 4x4 matrix from a linear std_vector_double
    mat=[]
    for iline in range(4):
        line=[]
        for icol in range(4):
            line.append(mat_std[iline*4+icol])
        mat.append(line)

    out=AttractRigidbody(rigid)
    for i in range(len(rigid)):
        coords=rigid.getCoords(i)
        coords2=Coord3D()
        coords2.x = mat[0][0]*coords.x + mat[0][1]*coords.y + mat[0][2]*coords.z + mat[0][3]
        coords2.y = mat[1][0]*coords.x + mat[1][1]*coords.y + mat[1][2]*coords.z + mat[1][3]
        coords2.z = mat[2][0]*coords.x + mat[2][1]*coords.y + mat[2][2]*coords.z + mat[2][3]
        out.SetCoords(i, coords2)
    return out

# check if a required file is found
def checkFile(name, comment):
    flag = os.path.exists(name)
    if not flag :
        msg =  "ERROR: file %s not found\n" %(name)
        if comment != "":
            msg += "ERROR: %s" %(comment)
        sys.exit(msg)	


###########################
##  MAIN ATTRACT PROGRAM  #
###########################
from optparse import OptionParser
parser = OptionParser(usage="%prog -r receptor_file -l ligand_file [-h] [-s] [-t] [--ref]")
parser.add_option("-r", "--receptor", action="store", type="string", dest="receptor_name", help="name of the receptor file")
parser.add_option("-l", "--ligand", action="store", type="string", dest="ligand_name", help="name of the ligand file")
parser.add_option("-s", "--single", action="store_true", dest="single", default=False, help="single minimization mode")
parser.add_option("--ref", action="store", type="string", dest="reffile", help="reference ligand for rmsd" )
parser.add_option("-t", "--translation", action="store", type="int", dest="transnb", help="translation number (distributed mode) starting from 0 for the first one!")
parser.add_option("--start1", action="store_true", default=False, dest="start1", help="(only useful with -t), use 1 for the first translation point")
(options, args) = parser.parse_args()


#receptor_name=args[0]
#ligand_name=args[1]

print """
**********************************************************************
**                                                                  **
**                ATTRACT  (Python edition)                         **
**                based on the PTools library                       **
**                                                                  **
**********************************************************************
PTools revision %s

"""%(Version().revid)

import locale

#locale.setlocale(locale.LC_ALL, 'fr_FR')
time_start = datetime.datetime.now()
#print now,"(",now.strftime("%A %B %d %Y, %H:%M"),")"
print "Start time:", time_start



#==========================
# read parameter file
#==========================

print "Reading parameters file: attract.inp"
(nbminim,lignames,minimlist,rstk) = readParams("attract.inp")
print "rstk = ",rstk


#open receptor and ligand files and check for the reduce format
def check_ffversion(reduced):
    header = open(reduced, 'r').readline()
    if not 'HEADER' in header:
         sys.stderr.write("ERROR: reduced PDB file must contain a HEADER line specifying the chosen forcefield (scorpion, attract1, attract2, imc)\n")
         sys.exit(1)

    #read cg format:
    return header.split()[1]
    
rec_ff = check_ffversion(options.receptor_name)
lig_ff = check_ffversion(options.ligand_name)

if rec_ff != lig_ff:
    sys.stderr.write("ERROR: reduction method differs between receptor and ligand\n")



allff_specs = {
             'SCORPION': {'ff_file': 'scorpion.par', 
                          'ff_class': ScorpionForceField,
                          'minimizer_class': ScorpionLbfgs
                          },

             'ATTRACT1': {'ff_file': 'aminon.par', 
                          'ff_class': AttractForceField1,
                          'minimizer_class': Lbfgs
                          },

             'ATTRACT2': {'ff_file': 'mbest1u.par', 
                          'ff_class': AttractForceField2,
                          'minimizer_class': Lbfgs
                          },
                          
             'IMC': {'ff_file': 'imc.par', 
                          'ff_class': ImcForceField,
                          'minimizer_class': ImcLbfgs
                          },
           }


ff_specs = allff_specs[rec_ff]


#==========================
# check required files
#==========================
# receptor
if not options.receptor_name:
    parser.print_help()
    parser.error("option -r is mandatory")
checkFile(options.receptor_name, "")
# ligand
if not options.ligand_name:
    parser.print_help()
    parser.error("option -l is mandatory")
checkFile(options.ligand_name, "")
# attract.inp
checkFile("attract.inp", "parameters file is required.")
# check if forcefield parameters file is present
checkFile(ff_specs['ff_file'], "forcefield file is required.")




#load receptor and ligand:
rec=Rigidbody(options.receptor_name)
lig=Rigidbody(options.ligand_name)
rec=AttractRigidbody(rec)
lig=AttractRigidbody(lig)
print "Reading receptor (fixed): %s with %d particules" %( options.receptor_name, len(rec) )
print "Reading  ligand (mobile): %s with %d particules" %( options.ligand_name,   len(lig) )

if (options.single and options.transnb):
    parser.error("options -s and -t are mutually exclusive")

# save all minimization variables in trajectory file
trjname = "minimization.trj"
if (options.single):
    ftraj = open(trjname, "w")

if (options.reffile):
    checkFile(options.reffile, "")
    ref=Rigidbody(options.reffile)
    print "Reading reference file: %s with %d particules" %( options.reffile, len(ref) )
    refca = ref.CA()
    if len(refca) == 0:  #No C alpha atom, ligand is probably a dna
        Rmsd_alias = Rmsd
        print "No Calpha atom found for ligand (DNA?). RMSD will be calculated on all grains"
    else:
        Rmsd_alias = rmsdca

if (not options.single):
    #systematic docking with default translations and rotations
    # check for rotation.dat and translation.dat
    checkFile("rotation.dat", "rotation file is required.")
    checkFile("translation.dat", "translation file is required.\nFormer users can rename translat.dat into translation.dat.")
    translations=Translation()
    rotations=Rotation()
else: #(single mode)
    #creates dummy translation and rotation
    translations=[[1,lig.FindCenter()]]
    rotations=[(0,0,0)]
    print "Single mode simulation"



printFiles=True
# option -t used: define the selected translation
transnb=0
if (options.transnb!=None):
    # check for rotation.dat and translation.dat
    checkFile("rotation.dat", "rotation file is required.")
    checkFile("translation.dat", "translation file is required.\nFormer users may rename translat.dat into translation.dat.")
    trans=Rigidbody("translation.dat")

    transnb=options.transnb

    if options.start1 is True:
       transnb -= 1

    co=trans.getCoords(transnb)
    translations=[[transnb+1,co]]

    if transnb!= len(trans)-1:
        printFiles=False #don't append ligand, receptor, etc. unless this is the last translation point of the simulation



# core attract algorithm
for trans in translations:
    transnb+=1
    print "@@@@@@@ Translation nb %i @@@@@@@" %(transnb)
    rotnb=0
    for rot in rotations:
        rotnb+=1
        print "----- Rotation nb %i -----"%rotnb
        minimcounter=0
        ligand=AttractRigidbody(lig)

        center=ligand.FindCenter()
        ligand.Translate(Coord3D()-center) #set ligand center of mass to 0,0,0
        ligand.AttractEulerRotate(surreal(rot[0]),surreal(rot[1]),surreal(rot[2]))
        ligand.Translate(trans[1])

        for minim in minimlist:
            minimcounter+=1
            cutoff=math.sqrt(minim['squarecutoff'])
            niter=minim['maxiter']
            print "{{ minimization nb %i of %i ; cutoff= %.2f (A) ; maxiter= %d"%(minimcounter,nbminim,cutoff,niter)


            #performs single minimization on receptor and ligand, given maxiter=niter and restraint constant rstk
            forcefield=ff_specs['ff_class'](ff_specs['ff_file'], surreal(cutoff)   )
            rec.setTranslation(False)
            rec.setRotation(False)
            
            forcefield.AddLigand(rec)
            forcefield.AddLigand(ligand)
            rstk=minim['rstk']  #restraint force
            #if rstk>0.0:
                #forcefield.SetRestraint(rstk)
            lbfgs_minimizer=ff_specs['minimizer_class'](forcefield)
            lbfgs_minimizer.minimize(niter)
            X=lbfgs_minimizer.GetMinimizedVars()  #optimized freedom variables after minimization


            #TODO: test and use CenterToOrigin() !
            output=AttractRigidbody(ligand)
            center=output.FindCenter()
            output.Translate(Coord3D()-center)
            output.AttractEulerRotate(surreal(X[0]), surreal(X[1]), surreal(X[2]))
            output.Translate(Coord3D(surreal(X[3]),surreal(X[4]),surreal(X[5])))
            output.Translate(center)

            ligand=AttractRigidbody(output)
            if (options.single):
                ntraj=lbfgs_minimizer.GetNumberIter()
                for iteration in range(ntraj):
                    traj = lbfgs_minimizer.GetMinimizedVarsAtIter(iteration)
                    for t in traj:
                        ftraj.write("%f "%t)
                    ftraj.write("\n")
                ftraj.write("~~~~~~~~~~~~~~\n")


        #computes RMSD if reference structure available
        if (options.reffile):
            rms=Rmsd_alias(ref, output)
        else:
            rms="XXXX"


        #calculates true energy, and rmsd if possible
        #with the new ligand position
        forcefield=ff_specs['ff_class'](ff_specs['ff_file'],  surreal(15))
        print "%4s %6s %6s %13s %13s"  %(" ","Trans", "Rot", "Ener", "RmsdCA_ref")
        pl = AttractPairList(rec, ligand,surreal(15))
        print "%-4s %6d %6d %13.7f %13s" %("==", transnb, rotnb, forcefield.nonbon8(rec,ligand,pl), str(rms))
        output.PrintMatrix()


#output compressed ligand and receptor:
if ( not options.single and printFiles==True): 
    print compress_file(options.receptor_name)
    print compress_file(options.ligand_name)
    print compress_file(ff_specs['ff_file'])
    print compress_file("translation.dat")
    print compress_file("rotation.dat")
    print compress_file("attract.inp")

# close trajectory file for single minimization 
if (options.single):
    ftraj.close()
    print "Saved all minimization variables (translations/rotations) in %s" %(trjname)

# print end and elapsed time
time_end = datetime.datetime.now()
#print "Finished at: ",now.strftime("%A %B %d %Y, %H:%M")
print "End time:", time_end
print "Elapsed time:", time_end - time_start

