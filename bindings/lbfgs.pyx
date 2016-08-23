from cython.operator cimport dereference as dereference
from libcpp.string cimport string

cdef extern from "lbfgs_interface.h" namespace "PTools":
    cdef cppclass CppLbfgs "PTools::Lbfgs":
        CppLbfgs(CppForceField&)
        void minimize(int)
        vector[double] GetMinimizedVars()
        vector[double] GetMinimizedVarsAtIter(int)
        int GetNumberIter()
        #void denormalize_weights()
        #void normalize_weights()
        vector[vector[double]] getWeights()

cdef class Lbfgs:
    
    cdef CppLbfgs * thisptr

    def __cinit__(self, forcefield):
        
        ff = <ForceField?> forcefield
        cdef CppForceField * ffptr = <CppForceField*?> ff.thisptr
        self.thisptr = new CppLbfgs(deref(ffptr))

    def __dealloc__(self):
        del self.thisptr

    def minimize(self, int maxiter):
        self.thisptr.minimize(maxiter)

    def GetMinimizedVars(self):
        cdef vector[double] vars = self.thisptr.GetMinimizedVars()
        out = []
        for i in xrange(vars.size()):
            out.append(vars[i])
        return out

    def GetNumberIter(self):
        return self.thisptr.GetNumberIter()

    def GetMinimizedVarsAtIter(self, int iter):
        cdef vector[double] vars = self.thisptr.GetMinimizedVarsAtIter(iter)
        out = []
        for i in xrange(vars.size()):
            out.append(vars[i])
        return out

    #def denormalize_weights(self):
    #    self.thisptr.denormalize_weights()

    #def normalize_weights(self):
    #    self.thisptr.normalize_weights()

    def getWeights(self):
        self.thisptr.getWeights()