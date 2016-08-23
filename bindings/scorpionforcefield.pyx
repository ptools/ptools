from cython.operator cimport dereference as deref
from libcpp.string cimport string


cdef extern from "forcefield.h" namespace "PTools":
    cdef cppclass CppForceField "PTools::ForceField":
        pass

cdef extern from "attractforcefield.h" namespace "PTools":
    cdef cppclass CppBaseAttractForceField "PTools::BaseAttractForceField":
        pass

cdef extern from "scorpionforcefield.h" namespace "PTools":
    cdef cppclass CppScorpionForceField "PTools::ScorpionForceField":
       CppScorpionForceField(string&, double)
       void AddLigand(CppAttractRigidbody&)
       double Function(vector[double]&)
       double getVdw()
       double getCoulomb()
       double nonbon8(CppAttractRigidbody& , CppAttractRigidbody& , CppAttractPairList & , int) 


cdef class ScorpionForceField (BaseAttractForceField):
   
    #cdef CppScorpionForceField* thisptr


    def __cinit__(self, filename, cutoff):
    
        # deallocate
        del self.thisptr
        self.thisptr = <CppForceField*> 0
        
        cdef char* c_filename
        cdef string * cppname

        c_filename = <char*> filename
        cppname = new string(c_filename)
        self.thisptr = <CppForceField *> new CppScorpionForceField(deref(cppname), cutoff)
        del cppname

    def __dealloc__(self):
        del self.thisptr

    def AddLigand(self, AttractRigidbody rig):
        cdef CppBaseAttractForceField* cpp_ptr = <CppBaseAttractForceField*> self.thisptr
        cpp_ptr.AddLigand(deref(<CppAttractRigidbody*>rig.thisptr))

    def Function(self, vec):
        cdef vector[double] v
        for el in vec:
           v.push_back(el)

        cdef CppBaseAttractForceField* cpp_ptr = <CppBaseAttractForceField*> self.thisptr
        return cpp_ptr.Function(v)
        
    def getVdw(self):
        cdef CppBaseAttractForceField* cpp_ptr = <CppBaseAttractForceField*> self.thisptr
        return cpp_ptr.getVdw()
    
    def getCoulomb(self):
        cdef CppBaseAttractForceField* cpp_ptr = <CppBaseAttractForceField*> self.thisptr
        return cpp_ptr.getCoulomb()

    def nonbon8(self, AttractRigidbody rec, AttractRigidbody lig, AttractPairList pl, verbose=False):
        cdef CppBaseAttractForceField* cpp_ptr = <CppBaseAttractForceField*> self.thisptr
        return cpp_ptr.nonbon8(deref(<CppAttractRigidbody*>rec.thisptr), deref(<CppAttractRigidbody*>lig.thisptr), deref(pl.thisptr), verbose)
