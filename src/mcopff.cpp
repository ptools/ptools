#include "mcopff.h"
#include <cassert>
#include <iostream>
#include <math.h>
#include <float.h>

namespace PTools
{

/////////////////// -- Class Mcoprigid -- ////////////////////

Mcoprigid::Mcoprigid(){
    _complete = false;

};

Mcoprigid::Mcoprigid(std::string filename){

    ReadMcoprigidPDB(filename);
    _center = _core.FindCenter();
    iniWeights();
    
}

// Normalized weight initialized at 1/(number of copies) and denormalized weights initialized at 0.
void Mcoprigid::iniWeights(){
    for(uint i=0; i < _vregion.size(); i++){
        std::vector<dbl> newvector;
        _weights.push_back(newvector);
        for(uint j=0; j < _vregion[i].size(); j++){
            dbl weight = 1.0/(double)_vregion[i].size();
            _weights[i].push_back(weight);
        }
    }
}

void Mcoprigid::setCore(AttractRigidbody& core) {
    _core=core;
    _center = _core.FindCenter();

};


void Mcoprigid::AttractEulerRotate(const dbl& phi, const dbl& ssi, const dbl& rot)
{
//Warning: makes euler rotation without centering
//the Mcoprigid object must be centered

//rotates the core body:
    _core.AttractEulerRotate(phi, ssi, rot);

//for each multicopy region, rotates the copy:
    for (uint i=0; i < _vregion.size(); i++)
        for (uint j=0; j<_vregion[i].size(); j++)
            _vregion[i][j].AttractEulerRotate(phi, ssi, rot);

}


void Mcoprigid::Translate(const Coord3D& c)
{

//translates the core body:
    _core.Translate(c);

//for each multicopy region, translates the copy:
    for (uint i=0; i < _vregion.size(); i++)
        for (uint j=0; j < _vregion[i].size(); j++)
            _vregion[i][j].Translate(c);

}


void Mcoprigid::ReadMcoprigidPDB(const std::string name) {

	// pointer toward the filename given in the constructor argument
    std::ifstream file(name.c_str()); 
    if (!file)
    {
        std::ostringstream oss;
        throw std::invalid_argument("##### ReadPDB:Could not open file \"" + name + "\" #####") ;
    }

    ReadMcoprigidPDB(file, _core, _vregion);
    file.close();

}


void Mcoprigid::ReadMcoprigidPDB(std::istream& file, AttractRigidbody& core, std::vector<AttractMcop>& regions){

    uint region_num = 0;
    uint copy_num = 0;
    std::string line;
    while(std::getline(file, line)){
        if(Mcop::isNewModel(line)){
            // The line is a new model
            if(line.size() > 12) region_num = line_to_region_number(line);
            if(line.size() > 15) copy_num = line_to_copy_number(line);
            AttractRigidbody model;
            while(std::getline(file,line)){
                if(isAtom(line)){
                    // The line is an atom
                    Coord3D pos = pdbToCoords(line);
                    Atomproperty a;
                    pdbToAtomproperty(line, a);
                    // if region_num is core, append attom to core
                    if(region_num == 0) core.AddAtom(a,pos);
                    else model.AddAtom(a,pos);
                }
                
                else{  
                    // if region_num is not core, add region copy
                    if(region_num != 0){
                        // if new region, add region to vector before adding region copy
                        if(regions.size() < region_num)
                            regions.push_back(AttractMcop());
                        // copy
                        regions[region_num-1].addCopy(model);
                    }
                    //Just finished adding a region copy
                    break;
                }
            }
        }
    }
}


void Mcoprigid::PrintWeights()
{
    for (uint loopregion=0; loopregion<_weights.size(); loopregion++)
    {
        std::cout << "copy region: " << loopregion << std::endl << "weights: ";
        for (uint copy=0; copy<_weights[loopregion].size(); copy++)
        {
            std::cout << _weights[loopregion][copy] << "  " ;
        }

    }
}

uint Mcoprigid::line_to_region_number(std::string line){

    return std::atoi(line.substr(12,1).c_str());    
}

uint Mcoprigid::line_to_copy_number(std::string line){

    return std::atoi(line.substr(15,line.size()-15).c_str());
}



/////////////////// -- Class Mcop -- ////////////////////

Mcop::Mcop(std::string filename){
    ReadmcopPDB(filename);
}


void Mcop::ReadmcopPDB(const std::string name) {

    std::string nomfich = name ;
	// pointer toward the filename given in the constructor argument
    std::ifstream file(nomfich.c_str()); 
    if (!file)
    {
        std::ostringstream oss;
        throw std::invalid_argument("##### ReadPDB:Could not open file \"" + nomfich + "\" #####") ;
    }

    ReadmcopPDB(file, _copies);
    file.close();

}

void Mcop::ReadmcopPDB(std::istream& file, std::vector<Rigidbody>& protein) {
    
    ReadModelsPDB(file, protein);

}

void Mcop::ReadModelsPDB(std::istream& file, std::vector<Rigidbody>& protein){

    std::string line;
    while(std::getline(file, line)){
        if(isNewModel(line)){
            // The line is a new model
            Rigidbody model;
            while(std::getline(file,line)){
                if(isAtom(line)){
                    // The line is an atom
                    Coord3D pos = pdbToCoords(line);
                    Atomproperty a;
                    pdbToAtomproperty(line, a);
                    model.AddAtom(a,pos);
                }
                else{
                    protein.push_back(model);
                    //Just finished adding a model
                    break;
                }
            }
        }
    }
}


bool Mcop::isNewModel(const std::string & line){

    if(line.substr(0,5)==(std::string)"MODEL") return true;
    else return false;
}


/////////////////// -- Class AttractMcop -- ////////////////////
AttractMcop::AttractMcop(std::string filename){
    
    Mcop copies(filename);
    for(uint i=0; i < copies.size(); i++){
        Rigidbody copy = copies.getCopy(i);
        AttractRigidbody attcopy = AttractRigidbody(copy);
        attract_copies.push_back(attcopy);
    }
    copies.clear();
}

AttractMcop::AttractMcop(const Mcop& mcop){
    
    Mcop copies = mcop;
    for(uint i=0; i < copies.size(); i++){
        Rigidbody copy = copies.getCopy(i);
        AttractRigidbody attcopy = AttractRigidbody(copy);
        attract_copies.push_back(attcopy);
    }
    copies.clear();
}


/////////////////// -- Class McopForceField -- ////////////////////

///////////////////////////////////////////////////
//     Forcefield implementation
///////////////////////////////////////////////////

void McopForceField::setReceptor( Mcoprigid& rec){

    Mcoprigid centeredrec = rec;
    //Coord3D center = rec.FindCenter();
    //rec._center = center;
    Coord3D center = _moved_ligand._center;

    centeredrec.Translate(Coord3D()-center);
    _receptor = centeredrec;
   
    ini_energies();
}

void McopForceField::setLigand( Mcoprigid& lig){

    Mcoprigid centeredlig = lig;
    Coord3D center = lig.FindCenter();
    lig._center = center;
    _moved_ligand = lig;

    centeredlig.Translate(Coord3D()-center);
    _centered_ligand = centeredlig;
}

void McopForceField::ini_energies(){
    for(uint i=0; i < _receptor.getRegions().size(); i++){
        std::vector<dbl> newvector;
        std::vector< std::vector<Coord3D> > newcoordsvecvec;
        std::vector<dbl> threezeros(3, 0.0);
        _Eik.push_back(newvector);
        _Z.push_back(0);
        _Zprime.push_back(0);
        _Eregion.push_back(0);
        _dEik.push_back(newcoordsvecvec);
        for(uint k=0; k < _receptor.getRegions()[i].size(); k++){
            _Eik[i].push_back(0);
            std::vector<Coord3D> newcoordsvec;
            _dEik[i].push_back(newcoordsvec);
        }
    }
}


uint McopForceField::ProblemSize()
{

    uint size = 0;
    if (_centered_ligand.getCore().hastranslation) size += 3;
    if (_centered_ligand.getCore().hasrotation) size += 3;
    return size;
}



/** \brief calculates energy of the system
*
* this functions returns nonbonded energy of a receptor with multicopy and a ligand without copy.
*
*/
dbl McopForceField::Function(const Vdouble & v)
{
    dbl ener_region = 0.0 ;
    dbl ener_core = 0.0 ;

// 1) put the objects to the right place
        
    _moved_ligand = _centered_ligand;

    Mcoprigid & lig = _moved_ligand ;
    assert(lig._vregion.size()==0);

    //TODO: take into account if no rotation or no translation
    lig.AttractEulerRotate(v[0],v[1],v[2]);
    lig.Translate(Coord3D(v[3],v[4],v[5]));
    assert(_receptor._vregion.size() == _receptor._weights.size());

//2) calculates the energy


    //2.1) core ligand body with core receptor
    AttractPairList pl (_receptor._core, lig._core, _cutoff );
    ener_core += _ff.nonbon8(_receptor._core, lig._core, pl );

    //2.2) core lignd with receptor copies:

    std::vector< std::vector<dbl> >& weights = _receptor._weights; 
    std::vector< std::vector<dbl> >& Eik = _Eik;
    std::vector< std::vector< std::vector<Coord3D> > >& dEik = _dEik;
    std::vector<dbl>& Z = _Z;
    std::vector<dbl>& Zprime = _Zprime;
    dbl& E = _E;
    std::vector<dbl>& Eregion = _Eregion;
    dbl& beta = _beta;

    for (uint loopregion=0; loopregion < _receptor._vregion.size() ; loopregion++)
    {

        //calculates interaction energy between receptor copies and ligand body:
        dbl enercopy =0.0;
        dbl& Z_i = Z[loopregion];
        dbl& Zprime_i = Zprime[loopregion];

        Z_i = 0.0;
        Zprime_i = 0.0;

        AttractMcop& ref_ensemble = _receptor._vregion[loopregion];
        std::vector<dbl>& weights_loop = weights[loopregion];
        
        assert( ref_ensemble.size() == weights_loop.size());

        double Emin = DBL_MAX;
        for (uint copynb = 0; copynb < ref_ensemble.size(); copynb++)
        {

            
            AttractRigidbody& copy = ref_ensemble[copynb];

            AttractPairList cpl (copy, lig._core, _cutoff );
            std::vector<Coord3D> copyforce(copy.Size());
            std::vector<Coord3D> coreforce(lig._core.Size());

            enercopy = _ff.nonbon8_forces(copy, lig._core, cpl, copyforce, coreforce);
            Eik[loopregion][copynb] = enercopy;
            if(enercopy < Emin) Emin = enercopy;
            
            //add force to ligand
            assert(lig._core.Size() == coreforce.size());
            
            dEik[loopregion][copynb] = coreforce;
        }

        for (uint copynb = 0; copynb < ref_ensemble.size(); copynb++){
            Z_i += exp(-beta*(Eik[loopregion][copynb] - Emin));
            Zprime_i += -Eik[loopregion][copynb]*exp(-beta*(Eik[loopregion][copynb] - Emin));
        }

        for (uint copynb = 0; copynb < ref_ensemble.size(); copynb++)
        {
            dbl& weight = weights_loop[copynb];
            weight = exp(-beta*(Eik[loopregion][copynb] - Emin))/Z_i;
        }
        Eregion[loopregion] = -Zprime_i/Z_i;
        ener_region += Eregion[loopregion];
    }
    E = ener_core + ener_region;
    return E;

}

// To calculation interaction energy between an Mcoprigid receptor with variable region(s) 
// and an Mcoprigid ligand without variable regions.
// Uses the receptor's normalized weights and a given forcefeild and cutoff. 

dbl McopForceField::CalcEnergy(Mcoprigid & receptor, Mcoprigid & lig, BaseAttractForceField & ff, dbl cutoff){

    dbl ener_region = 0.0 ;
    dbl ener_core = 0.0 ;

//2) calculates the energy

    //2.1) core ligand body with core receptor
    AttractPairList pl (receptor._core, lig._core, cutoff );
    ener_core += ff.nonbon8(receptor._core, lig._core, pl );

    //2.2) core lignd with receptor copies:

    std::vector< std::vector<dbl> >& weights = receptor._weights; 
    std::vector< std::vector<dbl> >& Eik = _Eik;
    std::vector<dbl>& Z = _Z;
    std::vector<dbl>& Zprime = _Zprime;
    dbl& E = _E;
    dbl beta = 1;

    for (uint loopregion=0; loopregion < receptor._vregion.size() ; loopregion++)
    {

        //calculates interaction energy between receptor copies and ligand body:

        dbl enercopy =0.0;
        dbl& Z_i = Z[loopregion];
        dbl& Zprime_i = Zprime[loopregion];

        Z_i = 0.0;
        Zprime_i = 0.0;


        AttractMcop& ref_ensemble = receptor._vregion[loopregion];
        std::vector<dbl>& weights_loop = weights[loopregion];
        
        assert( ref_ensemble.size() == weights_loop.size());


        for (uint copynb = 0; copynb < ref_ensemble.size(); copynb++)
        {

            AttractRigidbody& copy = ref_ensemble[copynb];

            AttractPairList cpl (copy, lig._core, cutoff );
            std::vector<Coord3D> copyforce(copy.Size());
            std::vector<Coord3D> coreforce(lig._core.Size());

            enercopy = ff.nonbon8_forces(copy, lig._core, cpl, copyforce, coreforce);
            Eik[loopregion][copynb] = enercopy;

            Z_i += exp(-beta*enercopy);
            Zprime_i += -enercopy*exp(-beta*enercopy);

        }
        for (uint copynb = 0; copynb < ref_ensemble.size(); copynb++)
        {
            dbl& weight = weights_loop[copynb];
            weight = exp(-beta*Eik[loopregion][copynb])/Z_i;
        }
        ener_region += -Zprime_i/Z_i;
    }
    E = ener_core + ener_region;
    return E;
}



void McopForceField::Derivatives(const Vdouble& v, Vdouble & g )
{

Mcoprigid & lig = _moved_ligand;
uint svptr = 0; // stateVars 'pointer'
// - force applied to the ligand by each loop copy
std::vector< std::vector< std::vector<Coord3D> > >& dEik = _dEik;
//Interaction nergy between each loop copy and the ligand
std::vector< std::vector<dbl> >& Eik = _Eik;
//Interactino energy between each variable region (sum of its copies) and the ligand.
std::vector<dbl>& Eregion = _Eregion;
//Weight of each loop copy
std::vector< std::vector<dbl> > & weights = _receptor._weights;
//beta coefficient (constant in energy function)
dbl& beta = _beta;

std::vector< std::vector<dbl> > coef;
std::vector<dbl> newvector;
// iterate through loop regions
for(uint i=0; i < dEik.size(); i++){
    coef.push_back(newvector);
    // iterate through loop copies
    for(uint k=0; k < dEik[i].size(); k++){
        coef[i].push_back((1+beta*(Eregion[i] - Eik[i][k]))*weights[i][k]);
    }
}

// calculate de rotational forces:
if (lig.getCore().hasrotation){
    //printf("hasrotation\n");
    dbl  cs,cp,ss,sp,cscp,sscp,sssp,crot,srot,xar,yar,cssp,x,y,z;
    dbl  pm[3][3];

    // !c     calculates orientational force contributions
    // !c     component 1: phi-angle
    // !c     component 2: ssi-angle
    // !c     component 3: rot-angle

    g[svptr+0] = 0.0;
    g[svptr+1] = 0.0;
    g[svptr+2] = 0.0;

    for (uint i=0; i<3;i++)
    {
        g[i+svptr]=0.0;
        for (uint j=0;j<3;j++)
            pm[i][j]=0.0 ;
    }

    cs=cos(v[svptr+1]);
    cp=cos(v[svptr+0]);
    ss=sin(v[svptr+1]);
    sp=sin(v[svptr+0]);
    cscp=cs*cp;
    cssp=cs*sp;
    sscp=ss*cp;
    sssp=ss*sp;
    crot=cos(v[svptr+2]);
    srot=sin(v[svptr+2]);

    //std::cout << "cs" << cs << std::endl;
    //std::cout << "cp" << cp << std::endl;
    //std::cout << "ss" << ss << std::endl;
    //std::cout << "sp" << sp << std::endl;
    //std::cout << "cscp" << cscp << std::endl;
    //std::cout << "cssp" << cssp << std::endl;
    //std::cout << "sscp" << sscp << std::endl;
    //std::cout << "sssp" << sssp << std::endl;
    //std::cout << "crot" << crot << std::endl;
    //std::cout << "srot" << srot << std::endl;

    // for the x, y and z coordinates, we need
    // the coordinates of the centered, non-translated molecule

    AttractRigidbody * pLigCentered = & _centered_ligand._core; // pointer to the centered ligand
    AttractRigidbody * pLigMoved  = & _moved_ligand._core; // pointer to the rotated/translated ligand (for forces)

    //initializing derivative variables
    std::vector<dbl> dEcore(3, 0.0);
    std::vector<dbl> dEloops(3, 0.0);

    for(uint atomIndex=0; atomIndex<_moved_ligand._core.Size(); atomIndex++){

        Coord3D coords = pLigCentered->GetCoords(atomIndex);

        x = coords.x;
        y = coords.y;
        z = coords.z;

        xar=x*crot+y*srot;
        yar=-x*srot+y*crot;

        pm[0][0]=-xar*cssp-yar*cp-z*sssp ;
        pm[1][0]=xar*cscp-yar*sp+z*sscp ;
        pm[2][0]=0.0 ;

        pm[0][1]=-xar*sscp+z*cscp ;
        pm[1][1]=-xar*sssp+z*cssp ;
        pm[2][1]=-xar*cs-z*ss ;

        pm[0][2]=yar*cscp+xar*sp ;
        pm[1][2]=yar*cssp-xar*cp ;
        pm[2][2]=-yar*ss ;

        //adding the derivatives for the core:
        for(uint j=0; j < 3; j++){
            dEcore[j] += pm[0][j] * pLigMoved->m_forces[atomIndex].x;
            dEcore[j] += pm[1][j] * pLigMoved->m_forces[atomIndex].y;
            dEcore[j] += pm[2][j] * pLigMoved->m_forces[atomIndex].z;
        }

        // iterate through loop regions
        for(uint i=0; i < dEik.size(); i++){
            // iterate through loop copies
            for(uint k=0; k < dEik[i].size(); k++){
                // iterate through the three rotation angles
                for(uint j=0; j < 3; j++){
                    dEloops[j] += coef[i][k] * pm[0][j] * dEik[i][k][atomIndex].x;
                    dEloops[j] += coef[i][k] * pm[1][j] * dEik[i][k][atomIndex].y;
                    dEloops[j] += coef[i][k] * pm[2][j] * dEik[i][k][atomIndex].z;
                }
            }
        }
    }


    for(uint j=0; j < 3; j++){
        g[svptr+j] += dEcore[j] + dEloops[j];
    }

    svptr += 3;


}

//sum the forces over x, y and z: calculate de translational forces
if (lig.getCore().hastranslation){

    g[svptr+0] = 0.0;
    g[svptr+1] = 0.0;
    g[svptr+2] = 0.0;

    dbl flim = 1.0e18;
    dbl ftr1, ftr2, ftr3, fbetr;
    ftr1=0.0;
    ftr2=0.0;
    ftr3=0.0;

    Coord3D ligtransForces; //translational forces for the ligand from loopregion
    
    //initializing derivative variables
    std::vector<dbl> dEcore(3, 0.0);
    std::vector<dbl> dEloops(3, 0.0);

    //adding derivatives
    for(uint atomIndex=0; atomIndex<_moved_ligand._core.Size(); atomIndex++){
        ligtransForces += _moved_ligand._core.m_forces[atomIndex];
        for(uint i=0; i < dEik.size(); i++){
            for(uint k=0; k < dEik[i].size(); k++){ 
                dEloops[0] += coef[i][k] * dEik[i][k][atomIndex].x;
                dEloops[1] += coef[i][k] * dEik[i][k][atomIndex].y;
                dEloops[2] += coef[i][k] * dEik[i][k][atomIndex].z;
            }
        }
    }

    dEcore[0] = ligtransForces.x;
    dEcore[1] = ligtransForces.y;
    dEcore[2] = ligtransForces.z;


    ftr1 = dEcore[0] + dEloops[0];
    ftr2 = dEcore[1] + dEloops[1];
    ftr3 = dEcore[2] + dEloops[2];
    //ftr1 =g[svptr+0];
    //ftr2 =g[svptr+1];
    //ftr3 =g[svptr+2];
   
    // force reduction, some times helps in case of very "bad" start structure
    for (uint i=0; i<3; i++)
    {
        fbetr=ftr1*ftr1 +ftr2*ftr2 +ftr3*ftr3;
        if (fbetr > flim)
        {
            ftr1=.01*ftr1;
            ftr2=.01*ftr2;
            ftr3=.01*ftr3;
        }
    }


    assert(svptr+2 < g.size());
    g[svptr+0] = ftr1;
    g[svptr+1] = ftr2;
    g[svptr+2] = ftr3;
    svptr += 3;
}


}


} // namespace PTools





