#include "lbfgs_interface.h"
#include "mcopff.h"

#include <string.h>
#include <iostream>
#include <math.h>
#include <float.h>

//#include "../complexify.h"

namespace PTools{

inline void assign(char* dest, char* src)
{
    memcpy(dest,src,strlen(src));

}






Lbfgs::Lbfgs( ForceField& toMinim)
        :objToMinimize(toMinim)
{
    //let the object do some initialization before beginning a new minimization
    //(for example, create new pairlists...)
    m_opt = NULL;
    objToMinimize.initMinimization();
};


Lbfgs::~Lbfgs()
{
    if (m_opt) lbfgsb_destroy(m_opt);
}


// vector<double> to vector<double> converter. used for genericity. should not impact performances too much.
inline void tocplx(const std::vector<double> & vdblin, std::vector<double> & vdblout ){vdblout=vdblin;};



#ifdef AUTO_DIFF
//convert from vector<cplx> to vector<double>
inline void tocplx(const std::vector<double> & vdblin, std::vector<surreal> & vcplx)
{
    vcplx = std::vector<surreal>();
    for (uint i=0; i<vdblin.size(); i++)
    {
        vcplx.push_back(vdblin[i]);
    }
}

#endif

inline std::vector<double> todbl(std::vector<double> & vdbl) {return vdbl;};

#ifdef AUTO_DIFF
inline std::vector<double> todbl(std::vector<surreal> & vcplx)
{
    std::vector<double> vdbl;
    for (uint i=0; i<vcplx.size(); i++)
    {
        vdbl.push_back(real(vcplx[i]));
    }
    return vdbl;
};

#endif

void Lbfgs::minimize(int maxiter)
{
    int n = objToMinimize.ProblemSize();
    std::cout  << "number of free variables for the minimizer: " << n << std::endl;
    std::vector<double> l(n);
    std::vector<double> u(n);
    Vint nbd(n);
    x.resize(n);
    g.resize(n);


    for (int i=0;i<n; i++)
    {
        l[i]=0;
        u[i]=0;
        nbd[i]=0;
        x[i] = 0.0;
        g[i] = 0.0;
    }  //unconstrained problem


    int rc;

    int m = 5;

    m_opt = lbfgsb_create(n, m, &l[0], &u[0], &nbd[0]);
    assert(m_opt);

    m_opt->iprint=-1;

    double f = DBL_MAX;

    m_opt->max_iter = maxiter;

    int last_iter = 0;
    int sub_iter = 0;

    double grad = 0;

    while (1) {
        rc = lbfgsb_run(m_opt, &x[0], &f, &g[0]);
        if (rc == 0) {
            break;
        } else if (rc < 0) {
            printf("lbfgsb stop with an error");
            break;
        } else if (rc == 1) {

            if (last_iter < m_opt->niter)
            {
                //this is a new iteration
                last_iter = m_opt->niter;
                sub_iter = 0;
                //saves the minimizer variables for each iteration (can be useful for generating animations)
                m_vars_over_time.push_back(x);
                objToMinimize.saveWeights(); //Only saves weights if objToMinimize is McopForceField. 
            }
            
            double sum = 0;
            for(int i=0; i < g.size(); i++){
                sum += g[i]*g[i];
            }
            

            std::vector<dbl> vdblx;
            tocplx(x,vdblx);
            std::vector<dbl> vdblg;
            tocplx(g,vdblg);
            f = objToMinimize.Function(vdblx);
            objToMinimize.Derivatives(vdblx,vdblg);
            g=todbl(vdblg);
            grad = sqrt(sum);
            sub_iter ++;

        } else {
            assert(!"can not reach here");
        }
    }


    m_opt->task[59]='\0'; //add a null terminating character to task[]


    std::cout << m_opt->task  << " |  " << m_opt->niter << " iterations\n";
}




std::vector<double>Lbfgs::GetMinimizedVarsAtIter(uint iter)
{
if (iter>=m_vars_over_time.size())
  {
   std::string msg = "";
   msg+= iter;
   msg += " is out of range (max: ";
   msg += m_vars_over_time.size()-1 ;
   msg += " )\n"; 
   throw std::out_of_range(msg);
  }
return m_vars_over_time[iter];
}

std::vector< std::vector<dbl> > Lbfgs::getWeights(){
    if (McopForceField * p = dynamic_cast<McopForceField *>(&objToMinimize)){
        // objToMinimize is of type McopForceField
        ForceField& r_objToMinimize = objToMinimize;
        McopForceField& r_Mcop_objToMinimize = dynamic_cast<McopForceField&>(r_objToMinimize);
        r_Mcop_objToMinimize.getWeights();
    }
}



} //namespace lbfgs






