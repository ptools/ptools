from cython.operator cimport dereference as deref
from libcpp.string cimport string


cdef extern from "atomselection.h" namespace "PTools":
    cdef cppclass CppAtomSelection "PTools::AtomSelection":
        CppAtomSelection()
        CppAtomSelection(CppAtomSelection&)
        CppAtomSelection(CppRigidbody&)

        unsigned int Size()
        void SetRigid(CppRigidbody&)
        CppAtom operator[]
        CppAtom CopyAtom(unsigned int)
        void AddAtomIndex(unsigned int)
        CppRigidbody CreateRigid()

        CppAtomSelection non(CppAtomSelection &)

    cdef CppAtomSelection operator& (CppAtomSelection&, CppAtomSelection&)
    cdef CppAtomSelection operator| (CppAtomSelection&, CppAtomSelection&)
    cdef CppAtomSelection operator! (CppAtomSelection&) 
#    cdef CppAtomSelection op_not (CppAtomSelection&)

cdef class AtomSelection:

    cdef CppAtomSelection * thisptr
    cdef object pyRigid

    def __cinit__(self, arg=None):
        cdef AtomSelection atsel
        cdef CppAtomSelection* atselptr
        cdef Rigidbody rig
        cdef CppRigidbody* rigptr

        self.thisptr = <CppAtomSelection*> 0
        if arg is None:
            self.thisptr = new CppAtomSelection()
            return

        elif isinstance(arg, AtomSelection):
            atsel  = <AtomSelection> arg
            atselptr  = <CppAtomSelection*> atsel.thisptr
            self.thisptr = new CppAtomSelection(deref(atselptr))
            return

        elif isinstance(arg, Rigidbody):
            rig = <Rigidbody> arg
            rigptr = rig.thisptr
            self.pyRigid = arg
            self.thisptr = new CppAtomSelection(deref(rigptr))

        else:
            raise RuntimeError("cannot reach here")

    def __dealloc__(self):
        if self.thisptr:
            del self.thisptr

    def __len__(self):
        return self.thisptr.Size()
    
    def __and__(AtomSelection self, AtomSelection second):
        ret = AtomSelection()
        del ret.thisptr
        cdef CppAtomSelection new_sel =   deref(self.thisptr) & deref(second.thisptr)    
        ret.thisptr  = new CppAtomSelection(new_sel)
        return ret

    def __or__(AtomSelection self, AtomSelection second):
        ret = AtomSelection()
        del ret.thisptr
        cdef CppAtomSelection new_sel =   deref(self.thisptr) | deref(second.thisptr)    
        ret.thisptr  = new CppAtomSelection(new_sel)
        return ret


    def __invert__(AtomSelection self):
        ret = AtomSelection()
        del ret.thisptr
        cdef CppAtomSelection new_sel =   not (deref(self.thisptr))  #uses C++ 'operator!'
        ret.thisptr  = new CppAtomSelection(new_sel)
        return ret

    def not_(self):
        return self.__invert__()
    
    
    def CreateRigid(self):
        ret = Rigidbody()
        if ret.thisptr:
            del ret.thisptr
        cdef CppRigidbody rig = self.thisptr.CreateRigid()
        ret.thisptr = new CppRigidbody(rig)
        return ret

    def SetRigid(self, Rigidbody r):
        self.pyRigid = r # to increase the refcount of r, preventing bad things if r is destroyed
        self.thisptr.SetRigid(deref(r.thisptr))
