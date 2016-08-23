from cython.operator cimport dereference as deref
from libcpp.string cimport string


cdef extern from "attractrigidbody.h" namespace "PTools":
    cdef cppclass CppAttractRigidbody "PTools::AttractRigidbody":
        CppAttractRigidbody()
        CppAttractRigidbody(string&)  #filename
        CppAttractRigidbody(CppRigidbody&)
        CppAttractRigidbody(CppAttractRigidbody&)

        unsigned int getAtomTypeNumber(unsigned int)

        double getCharge(unsigned int)

        int isAtomActive(unsigned int)

        void resetForces()
        #void addForces(  #TODO: later
 
        void setRotation(int)
        void setTranslation(int)
        CppCoord3D FindCenter()

        unsigned int Size()


        #returns radius of gyration
        double RadiusGyration()

        #returns the radius of a Rigidbody (max distance from center)
        double Radius()

        void PrintMatrix()

        CppAtomSelection CA()



cdef CppAttractRigidbody* _getAttractRigidbody_from_py_name(pyname):
    cdef char* name = pyname
    cdef string *cppname = new string(name)
    cdef CppAttractRigidbody *newrigid = new CppAttractRigidbody(deref(cppname))
    del cppname
    return newrigid


cdef class AttractRigidbody (Rigidbody) :

    def __cinit__(self, arg=''):
        
        # first deallocate the previously allocated Rigidbody
        del self.thisptr
        self.thisptr = <CppRigidbody*> 0

        if isinstance(arg, Rigidbody):
           rigidbody = <Rigidbody> arg
           rigidbodyptr = <CppRigidbody*> rigidbody.thisptr
           self.thisptr = <CppRigidbody*> new CppAttractRigidbody(deref(rigidbodyptr))
           return
        elif isinstance(arg, str):
           if arg == '':
              self.thisptr = <CppRigidbody*> new CppAttractRigidbody()
              return
           else:
              self.thisptr = <CppRigidbody*> _getAttractRigidbody_from_py_name(arg)
              return
        elif isinstance(arg, AttractRigidbody):
           oldrigidbody = <AttractRigidbody> arg
           oldrigidbody_ptr = <CppAttractRigidbody*> oldrigidbody.thisptr
           self.thisptr = <CppRigidbody*> new CppAttractRigidbody(deref(oldrigidbody_ptr) )          
           return
        else:
           ret =  "Should never reach here(attractrigidbody.pyx:AttractRigidbody:__cinit__)"
           print ret
           print arg
           raise ret


    def __dealloc__(self):
       if self.thisptr:
           del self.thisptr
           self.thisptr = <CppRigidbody*> 0

    
    def getAtomTypeNumber(self, atomid):
         return (<CppAttractRigidbody*>self.thisptr).getAtomTypeNumber(atomid)

    def getCharge(self, atomid):
         return (<CppAttractRigidbody*>self.thisptr).getCharge(atomid)

    #void setRotation(bool)
    def setRotation(self, flag):
        (<CppAttractRigidbody*> self.thisptr).setRotation(flag)

    #void setTranslation(bool)
    def setTranslation(self, flag):
        (<CppAttractRigidbody*> self.thisptr).setTranslation(flag)

    def isAtomActive(self, atomid):
        return (<CppAttractRigidbody*> self.thisptr).isAtomActive(atomid)

    #void resetForces()
    def resetForces(self):
        (<CppAttractRigidbody*> self.thisptr).resetForces()

    
    def Size(self):
        return self.thisptr.Size()

    #define also the __len__ method:
    def __len__(self):
        return self.thisptr.Size()

        
    def FindCenter(self):
        cdef CppRigidbody* rig = <CppRigidbody*> self.thisptr
        cdef CppCoord3D* co = new CppCoord3D (rig.FindCenter())
        ret = Coord3D()
        del ret.thisptr
        ret.thisptr = co
        
        return ret
        
    def Translate(self, Coord3D co):
        cdef CppRigidbody* rig = <CppRigidbody*> self.thisptr
        rig.Translate(deref(co.thisptr))
        
    def  AttractEulerRotate(self, double phi, double ssi, double rot):
        cdef CppRigidbody* rig = <CppRigidbody*> self.thisptr
        rig.AttractEulerRotate(phi, ssi, rot)


    #these function should be defined only in Rigdibody object and attractrigdbody should inherit from it:œ

    def Radius(self):
       return self.thisptr.Radius()

    def RadiusGyration(self):
       return self.thisptr.RadiusGyration()      


    def PrintMatrix(self):
       (<CppAttractRigidbody*> self.thisptr).PrintMatrix()

    def CA(self):
       ret = AtomSelection()
       del ret.thisptr
       cdef CppAtomSelection new_sel =  self.thisptr.CA()
       ret.thisptr  = new CppAtomSelection(new_sel)
       return ret
