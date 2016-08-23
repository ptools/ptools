#!/usr/bin/env python
# -*- coding: utf-8 -*-
from ptools import *
import sys


def multicopy(models_pdb, mcop_arg):
    """
    Takes the filename of a pdb containing multiple models and the string of
    the mcop argument indicating the region positions.

    Returns a single core of the protein (Rigidbody object) and all the copies
    of the region(s) (list of Mcop objects where each Mcop object contains the
    copies of a region which are Rigidbody) with respect to the mcop argument.

    ex: core, region_copies = ("NMRmodels.pdb", "23-44:66-79")

    The mcop argument uses the following format: 
        - If there are different regions, they are separated by ':'
        - The start and end position of a region are separated by '-'
        - To fetch copies of the whole protein, the mcop argument is ':'

    """
    starts, ends = convert_mcop_argument(mcop_arg)
    models = Mcop(models_pdb)
    cores, region_copies = extract_regions(models, starts, ends)
    if mcop_arg == ':':
        core = Rigidbody()
        moved_copies = region_copies
    else:
        core = find_core(cores)
        matrices = get_align_matrices(cores, core)
        moved_copies = []
        for i in xrange(0,len(region_copies)):
            moved_copies.append(apply_matrices(region_copies[i], matrices))
    return core, moved_copies


def convert_mcop_argument(arg):
    """
    Takes the mcop arguement (string).
    
    Returns start and end positions (list of int) of the region(s).

    ex: starts, ends = convert_mcop_argument("23-44:66-79")
    - starts variable contains [23, 66]
    - ends variable contains [66, 79]

    Note : if the argument is ':', returns (-1,) and (-1,). This case is for
    fetching copies of the whole protein. See extract_regions()
    """
    # if whole protein
    if arg == ':':
        return ((-1,), (-1,))
    else:
        reg = arg.split(':')
        starts = []
        ends = []
        for i in xrange(0,len(reg)):
            starts.append(min(map(int, reg[i].split('-'))))
            ends.append(max(map(int, reg[i].split('-'))))
        starts, ends = zip(*sorted(zip(starts,ends)))
        return starts, ends


def extract_regions(models, starts, ends):
    """
    Takes different models of a protein (Mcop object), start
    position(s) of region(s) (list int), and end position(s) of
    region(s) (list int).

    Returns the different cores of the protein (Mcop object)
    and all the copies of the region(s) (list of Mcop objects where each Mcop 
    object contains the copies of a region) with respect to the start and
    end positions.

    ex: cores, region_copies = extract_regions(models, starts, ends)

    Note: To fetch copies of the whole protein, starts and ends must both
    contain -1 (int) at index 0.
    
    Note: a single core can then be selected from the returned multiple cores
    thanks to the find_core() function.
    """
    mcop = []
    cores = Mcop()  
    # loop through models
    for i in xrange(0,len(models)):
        reg = []
        model = models[i]
        # selecting first region
        reg.append(model.SelectResRange(starts[0], ends[0]))
        fused_regions = reg[0]
        if i == 0:
            mcop.append(Mcop())
        #Creating and adding copy of region
        copy = reg[0]
        mcop[0].addCopy(copy.CreateRigid())
        for j in xrange(1, len(starts)):
            #adding other region
            if i == 0:
                mcop.append(Mcop())
            #Creating and adding copy of region
            reg.append(model.SelectResRange(starts[j], ends[j]))
            copy = reg[j]
            mcop[j].addCopy(copy.CreateRigid())
            fused_regions = copy | fused_regions
        #Creating and Adding core (inverse of regions)
        core = AtomSelection(model)
        core = core & ~fused_regions
        core = core.CreateRigid()
        cores.addCopy(core)
    # if whole protein
    if starts[0] == -1 and ends[0] == -1:
        temp = cores
        cores = mcop[0]
        mcop[0] = temp
    return cores, mcop


def write_region_copies(filename, models, append=False):
    """
    Writes a list of Mcop objects in a pdb file where each written model is
    delimited by MODEL i j and ENDMDL i j. (i = region number, j = copy number)

    Takes an output filename (string), and a list of Mcop objects.

    Default arguments : 
    - append option (False if overide, True if append), default value = False
    """
    if append:
        f = open(filename, 'a')
    else:
        f = open(filename, 'w')
    # iterate through number of regions
    for i in xrange(0, len(models)):
        # iterate through number of regions copies
        for j in xrange(0,len(models[i])):
            f.write('MODEL       ' + str(i+1) + '  ' + str(j+1) + '\n')
            f.write(str(models[i][j]) + '\n')
            f.write('ENDMDL    ' + str(i+1) + '  ' + str(j+1) + '\n')
    f.close()

    
def write_single_model(filename, model, append=False, model_index=1):
    """
    Writes a single Rigidbody object into a pdb file. The written model is
    delimited by MODEL i and ENDMDL i (i = index value of the model). 

    Takes an output filename (string), and a Rigidbody object.

    Default arguments : 
    - "append": False if overide file, True if append file, default value = False
    - "model_index": index value of the written model (int), default value = 1
    """
    if append:
        f = open(filename, 'a')
    else:
        f = open(filename, 'w')
    f.write('MODEL       ' + str(model_index) + '\n')
    f.write(str(model) + '\n')
    f.write('ENDMDL      ' + str(model_index) + '\n')
    f.close()

    
def write_mcop(filename, core, region_copies):
    """
    Writes a multicopy pdb file. This type of file uses the following format :
    The first model is the core of the protein and is delimited by MODEL 0 and
    ENDMDL 0.
    
    The following models are the copies of the defined variable region(s).
    These are delimited by MODEL i j and ENDMDL i j 
    (i = region number, j = copy number).

    Takes an output filename (string), the protein core (Rigidbody object),
    and the region copies (list of Mcop objects).
    """
    if len(core) > 0:
        write_single_model(filename, core, model_index=0)
    write_region_copies(filename, region_copies, append=True)


def rmsd_sums(cores):
    """
    Does pairwise superpositions of proteins, and for each protein, calculates
    the sum of the RMSD obtained from all of its superpositions. 

    Takes a list of proteins, generally protein cores 
    (list of Rigidbody objects).
    
    Returns a list of the corresponding RMSD sums (list of int).
    """
    rmsd_sums = [0 for x in xrange(0,len(cores))]
    for i in xrange(0, len(cores)):
        core1 = cores[i]
        for j in xrange(i+1,len(cores)):
            core2 = cores[j]
            sup = superpose(core1.CA().CreateRigid(), core2.CA().CreateRigid())
            rmsd_sums[i] += sup.rmsd
            rmsd_sums[j] += sup.rmsd
    return rmsd_sums


def find_core(cores):
    """
    From a list of proteins, returns the protein with the smallest RMSD sum
    (with pairwise superpositions using the rmsd_sums function). This
    is generally used to select the protein core from a list of potential cores. 

    Takes of a list of Rigidbody objects.
    
    Returns a single Rigidbody object.
    """
    sys.stderr.write("\nFinding core:\n")
    sums = rmsd_sums(cores)
    core_index = sums.index(min(sums))
    for i,sum in enumerate(sums): 
        sys.stderr.write("Core candidate %i: " %(i+1))
        sys.stderr.write("mean RMSD %f A \n" %(sum/(len(sums)-1)))
    sys.stderr.write("\nSelected core candiate %i: mean RMSD %f A\n"\
    %(core_index+1, sums[core_index]/(len(sums)-1)))
    return cores[core_index]


def get_align_matrices(cores, core):
    """
    Finds the displacement matrices (translation and rotation) of
    superpositions of a list of proteins (generally counterselected
    protein cores) over one protein (generally the selected protein core).
    These matrices are usefull to then move copies of other regions of the
    protein (they generally correspond to, but are not present in, the
    conterselected protein cores) to better position them onto the selected
    protein core. 

    Takes a Mcop objects, and a single Rigidbody objects. 
    Returns a list of Matrix objects corresponding to the entered list of
    Rigidbody objects. 
    """
    matrices = []
    for i in xrange(0, len(cores)):
        mobile_core = cores[i]
        sup = superpose(core, mobile_core)
        matrices.append(sup.matrix)
    return matrices


def apply_matrices(region_copies, matrices):
    """
    Modifies the coordinates of a list of Mcop objects according to a list of
    Matrix objects. This function is usually used after obtaining the matrices
    from the get_align_matrices function. 

    Takes a series of region copies (list of Mcop objects)
    Returns a series of  moved region copies (list of Mcop ojects with
    new cordinates)

    """
    moved_copies = Mcop()
    for i in xrange(0, len(region_copies)):
        cop = region_copies[i]
        cop.ApplyMatrix(matrices[i])
        moved_copies.addCopy(cop)
    return moved_copies

