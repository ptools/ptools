from cython.operator cimport dereference as deref
from libcpp.string cimport string


cdef extern from "imcforcefield.h" namespace "PTools":
    cdef cppclass CppImcForceField "PTools::ImcForceField":
       CppImcForceField(string&, double)
       void AddLigand(CppAttractRigidbody&)
       double Function(vector[double]&)
       double getVdw()
       double getCoulomb()
       double nonbon8(CppAttractRigidbody& , CppAttractRigidbody& , CppAttractPairList & , int) 


cdef class ImcForceField:
   
    cdef CppImcForceField* thisptr


    def __cinit__(self, filename, cutoff):
        cdef char* c_filename
        cdef string * cppname

        c_filename = <char*> filename
        cppname = new string(c_filename)
        self.thisptr = new CppImcForceField(deref(cppname), cutoff)
        del cppname

    def __dealloc__(self):
        del self.thisptr

    def AddLigand(self, AttractRigidbody rig):
        self.thisptr.AddLigand(deref(<CppAttractRigidbody*> rig.thisptr))

    def Function(self, vec):
        cdef vector[double] v
        for el in vec:
           v.push_back(el)

        return self.thisptr.Function(v)
        
    def getVdw(self):
        return self.thisptr.getVdw()
    
    def getCoulomb(self):
        return self.thisptr.getCoulomb()

    def nonbon8(self, AttractRigidbody rec, AttractRigidbody lig, AttractPairList pl, verbose=False):
        return self.thisptr.nonbon8(deref(<CppAttractRigidbody*>rec.thisptr), deref(<CppAttractRigidbody*>lig.thisptr), deref(pl.thisptr), verbose)
