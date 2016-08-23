#ifndef LBFGS_H
#define LBFGS_H

#include "basetypes.h"
#include "forcefield.h"


#include "lbfgsb.h"




namespace PTools
{



// new version, from sdrive.c
// no bounds !!
class Lbfgs
{
      public:
            Lbfgs(ForceField& toMinim);
            ~Lbfgs();
            void minimize(int maxiter);
            std::vector<double> GetMinimizedVars() const {return x;};

            std::vector<double> GetMinimizedVarsAtIter(uint iter);
            int GetNumberIter() {return m_opt->niter;}
            //void denormalize_weights();
            //void normalize_weights();
            std::vector< std::vector<dbl> > getWeights();




      private:

            ForceField& objToMinimize ;
            std::vector<double> x ; // position variables
            std::vector<double> g ; // gradient

            lbfgsb_t* m_opt; //minimizer structure

            std::vector<std::vector<double> > m_vars_over_time;



} ;

}

#endif //#ifndef Lbfgs_H

