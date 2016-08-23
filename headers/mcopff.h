//  multicopy-forcefield
//
//



#ifndef _MCOPFF_H_
#define _MCOPFF_H_

#include "attractforcefield.h"
#include "rigidbody.h"
#include "pdbio.h" 

namespace PTools{


// typedef std::vector<AttractRigidbody> ensemble ;

class Mcop
{
private:

std::vector<Rigidbody> _copies;

public:

     //Contructor : no arguments
     Mcop(){ }
     //Contructor : takes filename
     Mcop(std::string filename);
     
     
    virtual void addCopy(const Rigidbody& cop){_copies.push_back(cop);};

    
    virtual size_t size() const {return _copies.size();};
     

    virtual Rigidbody& operator[](uint i){return _copies[i];};
    virtual Rigidbody& getCopy(int i){return _copies[i];};

    void ReadmcopPDB(const std::string name);
    void ReadmcopPDB(std::istream& file, std::vector<Rigidbody>& protein);
    void ReadModelsPDB(std::istream& file, std::vector<Rigidbody>& protein);
    static bool isNewModel(const std::string & line);  
    
    void clear(){_copies.clear();};

};



class AttractMcop: public Mcop
{
private:

    std::vector<AttractRigidbody> attract_copies;

public:

    //Constructor: no arguments
    AttractMcop(){ };
    //Constructor: takes filename
    AttractMcop(std::string filename);
    //Constructor: takes mcop
    AttractMcop(const Mcop& mcop);
    
    virtual void addCopy(const AttractRigidbody& cop){attract_copies.push_back(cop);};

    virtual size_t size() const {return attract_copies.size();};
     

    virtual AttractRigidbody& operator[](uint i){return attract_copies[i];};
    virtual AttractRigidbody& getCopy(int i){return attract_copies[i];};
    
};



class Mcoprigid //multicopy rigidbody
{

public:
    //Constructor: no arguments
    Mcoprigid();
    //Constructor: takes filename
    Mcoprigid(std::string filename);
    //using default copy operator

    void iniWeights();
    void setCore(AttractRigidbody& core) ;
    void addEnsemble(const AttractMcop& reg){ _vregion.push_back(reg); std::vector<dbl> v; _weights.push_back(v);  };


    void AttractEulerRotate(const dbl& phi, const dbl& ssi, const dbl& rot);
    void Translate(const Coord3D& c);
    Coord3D FindCenter(){return _core.FindCenter();};

    void PrintWeights();
    std::vector <std::vector<dbl> > getWeights(){return _weights;};
    
    void ReadMcoprigidPDB(const std::string name);
    void ReadMcoprigidPDB(std::istream& file, AttractRigidbody& core, std::vector<AttractMcop>& regions);
    uint line_to_region_number(std::string line);
    uint line_to_copy_number(std::string line);
    
    AttractRigidbody& getCore(){return _core;};
    std::vector<AttractMcop>& getRegions(){return _vregion;};
    AttractMcop& getRegion(int i){return _vregion[i];};
    size_t size() const {return _vregion.size();};

    virtual void setRotation(bool value) {_core.setRotation(value);}; ///< allow/disallow rotation
    virtual void setTranslation(bool value) {_core.setTranslation(value);}; ///< allow/disallow translation
    bool checkRotation(){return _core.checkRotation();};
    bool checkTranslation(){return _core.checkTranslation();};

private:

    AttractRigidbody _core;
    std::vector< AttractMcop > _vregion ;

    bool _complete ; //is this used ?
    Coord3D _center ; ///<center of mass of the core region
    std::vector< std::vector<dbl> > _weights;

    friend class McopForceField;

};



/** \brief ForceField with multicopy
needs an attract forcefield (either 1 or 2) in constructor

no normal modes yet
*/
class McopForceField: public ForceField
{

public:

    McopForceField(BaseAttractForceField& ff, dbl cutoff)
            :_ff(ff), _cutoff(cutoff) { _beta = 0.3; };
    void ini_energies();

    dbl Function(const Vdouble&);
    dbl CalcEnergy(Mcoprigid & receptor,Mcoprigid & lig, BaseAttractForceField & ff, dbl cutoff);
    void Derivatives(const Vdouble& v, Vdouble & g );

    void setReceptor(Mcoprigid& rec);
    void setLigand(Mcoprigid& lig);

    uint ProblemSize();
    void initMinimization(){};
    virtual void saveWeights(){_savedWeights.push_back(_receptor.getWeights());};
    std::vector< std::vector< std::vector<dbl> > > getSavedWeights(){return _savedWeights;};

    std::vector <std::vector<dbl> > getWeights(){return _receptor.getWeights();};
    std::vector <std::vector<dbl> > getMcopE(){return _Eik;};


private:

    BaseAttractForceField& _ff ;
    dbl _cutoff;

    bool _update_weights;

    Mcoprigid _centered_ligand ;
    Mcoprigid _moved_ligand ;
    Mcoprigid _receptor;
    std::vector< std::vector<dbl> > _Eik; //Loop copy interaction energies
    std::vector< std::vector< std::vector<Coord3D> > > _dEik; //- force applied to the ligand by each loop copy
    dbl _E; //Total interaction energy
    std::vector<dbl> _Eregion; //Loop region interaction energies
    std::vector<dbl> _Z; // one Z per loop region
    std::vector<dbl> _Zprime; // one Z' per loop region
    dbl _beta; // constant in the energy function 
    std::vector< std::vector< std::vector<dbl> > > _savedWeights; //conserved weights throughout minimizations of a docking


    


};

}//namespace PTools

#endif // _MCOPFF_H_

