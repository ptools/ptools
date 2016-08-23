#ifndef FORCEFIELD_H
#define FORCEFIELD_H


#include "basetypes.h"
#include "rigidbody.h"
#include "attractrigidbody.h"
#include "pairlist.h"

namespace PTools
{


void PrintVec(const Vdouble& vec);


//! Generic forcefield (abstract class)
/**

 */
class ForceField
{
public:

    ///the function to optimize
    virtual dbl Function(const Vdouble&)=0;

    /// analytical derivative of the function above
    virtual void Derivatives(const Vdouble& StateVars, Vdouble& delta)
    {
        NumDerivatives(StateVars, delta, true);
    }

    ///numerical derivative for testing purpose. Not very accurate
    virtual void NumDerivatives(const Vdouble& StateVars, Vdouble& delta, bool print=false);

    ///size of the problem (number of variables the minimizer must optimize)
    virtual uint ProblemSize()=0;

    ///this function is called at the beginning of a minimization, by the minimizer (Lbfgs)
    virtual void initMinimization()=0;

    ///virtual destructor (Effective C++ - Scott Meyers - Item 7)
    virtual ~ForceField(){};
    
    /// Saves weights into attribute. Used at each minimization iteration in order to follow
    /// the weights' evolution during the docking. Used only for McopForcefield. 
    virtual void saveWeights(){};

} ;



} //namespace PTools





#endif //#ifndef FORCEFIELD_H

