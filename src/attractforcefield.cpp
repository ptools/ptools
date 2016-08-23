#include "attractforcefield.h"


#include <fstream>
#include <math.h>  //for fabs()
#include <sstream> //for istringstream


using std::ios_base;


namespace PTools {

/*! \brief Extracts extra information from ATOM lines.
 *
 *  For Attract pdb files, the library currently reads some extra informations
*  after the x,y,z coordinates. This information is extracted here.
*  Two arrays are populated then: atcategory which contains the atom type category
* (AKA iaci variable in the fortran code) and the atom charge ('xlai' in the fortran code)
 */
void extractExtra( Rigidbody& rig, std::vector<uint>& vCat, std::vector<dbl>& vCh)
{

    uint   atcategory  = 0;
    dbl  atcharge   = 0.0;

    for (uint i=0; i<rig.Size(); i++)
    {
        const Atomproperty at = rig.GetAtomProperty(i);
        std::string extra = at.extra;
        //std::cout << extra << std::endl;
        std::istringstream iss( extra );
        iss >> atcategory >> atcharge ;
        vCat.push_back(atcategory-1);
        vCh.push_back(atcharge);
    }


}



AttractForceField1::AttractForceField1(std::string paramsFileName, dbl cutoff)
//         :m_refreceptor(recept), m_refligand(lig), m_receptor(recept), m_ligand(lig),m_savligand(lig),
//         plist(recept,lig,cutoff)

{

    InitParams(paramsFileName);

    m_rstk=0.0; //no restraint by default
    m_cutoff = cutoff;

}



void AttractForceField1::InitParams(const std::string & paramsFileName )
{
    int indice, inull;
    dbl rad;
    dbl amp;

    std::ifstream aminon(paramsFileName.c_str());
    if (!aminon)
    {
        //the file cannot be opened
        std::string msg = "Forcefield.cpp: Cannot Locate file forcefield parameters (aminon.par)\n";
        std::cout << msg ;
        throw std::invalid_argument(msg);
    }

    while (!aminon.eof())
    {
        aminon >> indice >> rad >> amp >> inull ;
        m_rad.push_back(rad) ;
        m_amp.push_back(amp) ;
        assert(m_rad.size()<64);
    }


    //initialisation of the pre-calculated array of rc and ac
    for (uint i=0; i<m_rad.size();i++)
        for (uint j=0; j<m_rad.size(); j++)
        {
            m_rc[i][j]=m_amp[i]*m_amp[j]*pow((m_rad[i]+m_rad[j]),8);
            m_ac[i][j]=m_amp[i]*m_amp[j]*pow((m_rad[i]+m_rad[j]),6);
        }
}



dbl AttractForceField1::nonbon8_forces(AttractRigidbody& rec, AttractRigidbody& lig, AttractPairList & pairlist, std::vector<Coord3D>& forcerec, std::vector<Coord3D>& forcelig, bool print)
{

    assert(forcerec.size() == rec.Size());
    assert(forcelig.size() == lig.Size());

    dbl sumLJ=0.0 ;
    dbl sumElectrostatic=0.0;


    //synchronize coordinates for using unsafeGetCoords
    rec.syncCoords();
    lig.syncCoords();

    Coord3D a, b;


    for (uint iter=0; iter<pairlist.Size(); iter++)
    {

        uint ir = pairlist[iter].atrec;
        uint jl = pairlist[iter].atlig;

        uint rAtomCat = rec.getAtomTypeNumber(ir);
        uint lAtomCat = lig.getAtomTypeNumber(jl);

        assert(rAtomCat < m_rad.size());
        assert(lAtomCat < m_rad.size());


        dbl alen = m_ac[ rAtomCat ][ lAtomCat ];
        dbl rlen = m_rc[ rAtomCat ][ lAtomCat ];


        lig.unsafeGetCoords(jl,a);
        rec.unsafeGetCoords(ir,b);

        Coord3D dx = a-b ;
        dbl r2 = Norm2(dx);

        if (r2 < 0.001 ) r2=0.001;
        dbl rr2 = 1.0/r2;
        dx = rr2*dx;

        dbl rr23 = rr2*rr2*rr2 ;
        dbl rep =  rlen*rr2 ;
        dbl vlj = (rep-alen)*rr23 ;

        sumLJ += vlj;

        dbl fb = 6.0*vlj+2.0*(rep*rr23) ;
        Coord3D fdb = fb*dx ;

        //assign force to the atoms:
        forcelig[jl] -= fdb ;
        forcerec[ir] += fdb ;



        //electrostatic part:
        dbl chargeR = rec.m_charge[ir];
        dbl chargeL = lig.m_charge[jl];
        dbl charge = chargeR * chargeL * (332.053986/20.0);

        if (fabs(charge) > 0.0)
        {

            dbl et = charge*rr2;
            sumElectrostatic+=et;

            Coord3D fdb = (2.0*et)*dx;
            forcelig[jl] -= fdb ;
            forcerec[ir] += fdb ;
        }
    }

    m_vdw = sumLJ;
    m_elec = sumElectrostatic;

    return sumLJ + sumElectrostatic;
}









void BaseAttractForceField::initMinimization()
{
    MakePairLists();
}


////////////////////////////////////////////////////////////////
//     AttractForceField2 implementation
////////////////////////////////////////////////////////////////


static AttFF2_params* m_params = 0;

AttractForceField2::AttractForceField2(const std::string & filename, dbl cutoff)
{
     loadParams(filename, cutoff);
}

void AttractForceField2::resetParams()
{
     delete m_params;
     m_params = NULL;
     m_filename = "";
}

void AttractForceField2::reloadParams(const std::string & filename, dbl cutoff)
{
  resetParams();
  loadParams(filename, cutoff);
}


void AttractForceField2::loadParams(const std::string & filename, dbl cutoff)
{

    m_cutoff=cutoff;
    if (m_params==0)
    {
        m_params=new AttFF2_params();

        std::ifstream mbest (filename.c_str());
        //open(11,file=eingabe2) -> eingabe2: mbest1k.par
        if (!mbest)
        {
            //the file cannot be opened
            std::string msg = "Forcefield.cpp: Cannot Locate file  " + filename + "\n" ;
            ios_base::failure fail(msg);
            std::cout << msg ;
            throw fail;
        }
        m_filename = filename;



        std::string line;
        getline(mbest, line); //read the first line into "line"
        std::istringstream iss (line); // iss allows formated extraction from "line"

        std::string magic;
        iss >> magic;   //reads magic constant
        if (magic == "AFF")   // file format version 2 at least
        {
           int revnb;
           iss >> revnb;  //get revision number
           if (revnb == 2 ) // file format version 2:
             {
                 //get list of dummy types:
                 uint numdummy;
                 iss >> numdummy;
                 std::vector<uint> dummyatomtypes;
                 for(uint i=0; i<numdummy; i++)
                    {
                       uint type;
                       iss >> type;
                       dummyatomtypes.push_back(type-1); //types counting begins at 0
                    }
                 std::swap(dummyatomtypes, m_params->_dummytypes);
             }
            else
             {
                std::string msg = "AttractForceField2: cannot read parameter file version \n";
                msg += revnb ;
                std::cerr << msg ;
                ios_base::failure fail(msg);
                throw fail;
             }
        }
        else 
        {
          //declare forcefield file as invalid:
          std::string msg = "AttractForceField2: invalid paramters file format: doesn't contain AFF string \n";
          std::cerr << msg ;
          ios_base::failure fail(msg);
          throw fail;
        }


        for (uint i = 0; i<31; i++)
            for (uint j = 0; j<31; j++)
            {
                mbest >> m_params->rbc[i][j] ;
            }

        for (uint i = 0; i<31; i++)
            for (uint j = 0; j<31; j++)
                mbest >> m_params->abc[i][j] ;

        for (uint i = 0; i<31; i++)
        {
            for (uint j = 0; j<31; j++)
            {
                mbest >> m_params->iflo[i][j] ;
                assert(m_params->iflo[i][j]==1 || m_params->iflo[i][j]==-1);
            }
        }



        for (uint jj=0; jj<31; jj++)  // loop over attract atom types
        {

            for (uint ii=0; ii<31; ii++) // loop over attract atom types
            {

                dbl rbc2 = m_params->rbc[ii][jj]*m_params->rbc[ii][jj];
                dbl rbc6 = rbc2*rbc2*rbc2;
                dbl rbc8 = rbc6*rbc2;
                m_params->rc[ii][jj] = m_params->abc[ii][jj] * rbc8;
                m_params->ac[ii][jj] = m_params->abc[ii][jj] * rbc6;

                m_params->ipon[ii][jj] = m_params->iflo[ii][jj] ;
                assert(m_params->ipon[ii][jj]==1 || m_params->ipon[ii][jj]==-1);

                dbl alen = m_params->ac[ii][jj];
                dbl rlen = m_params->rc[ii][jj];
                dbl alen4 = alen*alen*alen*alen;
                dbl rlen3 = rlen*rlen*rlen;
                m_params->emin[ii][jj] = -27.0*alen4/(256.0*rlen3);
                m_params->rmin2[ii][jj]= 4.0*rlen/(3.0*alen);


            }
        }
    }


}


dbl BaseAttractForceField::Function(const Vdouble& stateVars )
{

    assert(m_centeredligand.size() >=1);
    assert(m_movedligand.size() >=1);
    
    if (stateVars.size() != this->ProblemSize() )
    {
      throw std::runtime_error("error: ProblemSize != size of state vars in BaseAttractForceField::Function");
    }

    uint svptr = 0; //state variable 'pointer'
    const uint nlig = m_movedligand.size();


    //don't let the user call this function without a coherent pairlist
    //(the pairlist may be outdated (user choice), but we MUST have the correct number of pairlists!)
    if (m_pairlists.size() != (nlig*(nlig-1))/2)
        MakePairLists();

    assert(m_pairlists.size() == (nlig*(nlig-1))/2);


    //put the ligands to the correct positions defined by stateVars
    for (uint i=0; i<m_movedligand.size(); i++)
    {

        m_movedligand[i] = m_centeredligand[i];
        m_movedligand[i].resetForces(); //just to be sure that the forces are set to zero. Maybe not needed.

        

        if (m_movedligand[i].hasrotation)
        {
            assert(svptr+2 < stateVars.size());
            m_movedligand[i].AttractEulerRotate(stateVars[svptr], stateVars[svptr+1], stateVars[svptr+2]);
            svptr+=3;
        }


        m_movedligand[i].Translate(m_ligcenter[i]);

        if (m_movedligand[i].hastranslation)
        {
            assert(svptr+2 < stateVars.size());
            m_movedligand[i].Translate(Coord3D(stateVars[svptr],stateVars[svptr+1],stateVars[svptr+2]));
            svptr+=3;
        }



    }




    dbl enernon = 0.0 ;

    uint plistnumber = 0; //index of pairlist used for a given pair of ligands
    //iteration over all ligand pairs:
    for (uint i=0; i<m_movedligand.size(); i++)
        for (uint j=i+1; j<m_movedligand.size(); j++)
        {
            assert(plistnumber < m_pairlists.size() );
            enernon += nonbon8(m_movedligand[i], m_movedligand[j],  m_pairlists[plistnumber++] );   //calculates energy contribution for every pair. Forces are stored for each ligand
        }


    return enernon;

}



uint BaseAttractForceField::ProblemSize()
{
    uint size = 0;
    for (uint i = 0; i < m_centeredligand.size(); i++)
    {
        if (m_centeredligand[i].hastranslation) size +=3 ;
        if (m_centeredligand[i].hasrotation) size +=3 ;
    }

    return size;
}






/*! \brief returns the analytical derivatives of the forcefield 2
*
*   input:
*   Vdouble & stateVars: determines how the molecules are moved by the minimizer
*   (the minimizer only works on a linear Vdouble holding the free minimization variables)
*   output:
*   this function puts the derivative of the energy with respect to variable 1 to 6 (in case of 6 degrees of
*   freedom. 3 trans + 3 rotations) into the 'delta' Vdouble array
*/
void BaseAttractForceField::Derivatives(const Vdouble& stateVars, Vdouble& delta)
{

    uint svptr = 0; // stateVars 'pointer'


    for (uint i=0; i<m_movedligand.size(); i++)
    {

        if (m_movedligand[i].hasrotation)
        {
            //calculates the rotational force for ligand i
            Rota(i, stateVars[svptr], stateVars[svptr+1], stateVars[svptr+2], delta, svptr, false );
            svptr+=3;
        }

        if (m_movedligand[i].hastranslation)
        {
            //calculates the translational force for ligand i
            Trans(i, delta, svptr, false);
            svptr+=3;
        }


    }

}






/*! \brief Non bonded energy
*
*   translated from fortran file nonbon8.f
*   TODO: add comments in the code, remove debug instructions
*/
dbl AttractForceField2::nonbon8_forces(AttractRigidbody& rec, AttractRigidbody& lig, AttractPairList & pairlist, std::vector<Coord3D>& forcerec, std::vector<Coord3D>& forcelig, bool print)
{

    dbl enon = 0.0;
    dbl epote = 0.0;


    std::cout.precision(20);

    //synchronise coordinates to later use unsafeGetCoords (should be faster)
    rec.syncCoords();
    lig.syncCoords();

    Coord3D a;
    Coord3D b;

    for (uint ik=0; ik<pairlist.Size(); ik++ )
    {
        AtomPair atpair = pairlist[ik];


        uint i = atpair.atrec ;
        uint j = atpair.atlig ;
        uint ii=rec.m_atomTypeNumber[i];
        uint jj=lig.m_atomTypeNumber[j];

        assert(ii<31);
        assert(jj<31);
        dbl alen = m_params->ac[ii][jj];
        dbl rlen = m_params->rc[ii][jj];
        int ivor = m_params->ipon[ii][jj];
        assert(ivor==1 || ivor==-1);


        dbl charge= rec.m_charge[i]* lig.m_charge[j];  //charge product of the two atoms

        rec.unsafeGetCoords(i,a); lig.unsafeGetCoords(j,b);

        Coord3D dx  ( a-b ) ;


        dbl r2 = Norm2(dx);
        if (r2 < 0.001) r2=0.001 ;

        dbl rr2 = 1.0/r2;
        dx = rr2*dx ;

        if (charge != 0.0) {
            dbl et = charge*rr2;
            et*=(332.053986/15.0);  //constant felec/permi (could still be optimized!)

            epote += et ;
            Coord3D fdb =2.0*et*dx ;
            forcelig[j] += fdb;
            forcerec[i] -= fdb;

        }

        //switch between minimum or saddle point
        if (r2 < m_params->rmin2[ii][jj] ) {

            dbl rr23 = rr2*rr2*rr2 ;
            dbl rep = rlen*rr2 ;
            dbl vlj = (rep-alen)*rr23;
            enon=enon+vlj+(ivor-1)*m_params->emin[ii][jj] ;

            dbl fb=6.0*vlj+2.0*(rep*rr23);
            Coord3D fdb = fb*dx;
            forcelig[j]+=fdb;
            forcerec[i]-=fdb;
        }
        else {
            dbl rr23=rr2*rr2*rr2;
            dbl rep=rlen*rr2;
            dbl vlj=(rep-alen)*rr23 ;
            enon += ivor*vlj ;

            dbl fb=6.0*vlj+2.0*(rep*rr23);

            Coord3D fdb=ivor*fb*dx ;
            forcelig[j] += fdb ;
            forcerec[i] -= fdb ;
        }


    }

    if (print) std::cout << "vlj  coulomb: " << enon << "  " << epote << "\n";
    m_elec = epote;
    m_vdw = enon;
    return enon+epote;
}



void BaseAttractForceField::Trans(uint molIndex, Vdouble & delta, uint shift,  bool print)
{
// molIndex is the index of the protein we want to extract the average
// translational forces


    AttractRigidbody const & rig(m_movedligand[molIndex]);
//   In this subroutine the translational force components are calculated
    dbl flim = 1.0e18;
    dbl ftr1, ftr2, ftr3, fbetr;

    ftr1=0.0;
    ftr2=0.0;
    ftr3=0.0;
    for (uint i=0;i<rig.Size(); i++)
    {
        ftr1=ftr1 + rig.m_forces[i].x;
        ftr2=ftr2 + rig.m_forces[i].y;
        ftr3=ftr3 + rig.m_forces[i].z;
    }

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

    assert(shift+2 < delta.size());
    delta[0+shift]=ftr1;
    delta[1+shift]=ftr2;
    delta[2+shift]=ftr3;

    //debug:
    if (print) std::cout <<  "translational forces: " << ftr1 <<"  "<< ftr2 <<"  " << ftr3 << std::endl;
    return ;
}



void BaseAttractForceField::Rota(uint molIndex, dbl phi,dbl ssi, dbl rot, Vdouble & delta,uint shift, bool print)
{
// molIndex is the index of the protein we want to extract the average
// translational forces



    //delta array of dbls of dimension 6 ( 3 rotations, 3 translations)

    dbl  cs,cp,ss,sp,cscp,sscp,sssp,crot,srot,xar,yar,cssp,X,Y,Z ;
    dbl  pm[3][3];

// !c
// !c     calculates orientational force contributions
// !c     component 1: phi-angle
// !c     component 2: ssi-angle
// !c     component 3: rot-angle
// !c


    for (uint i=0; i<3;i++)
    {
        delta[i+shift]=0.0;
        for (uint j=0;j<3;j++)
            pm[i][j]=0.0 ;
    }

    cs=cos(ssi);
    cp=cos(phi);
    ss=sin(ssi);
    sp=sin(phi);
    cscp=cs*cp;
    cssp=cs*sp;
    sscp=ss*cp;
    sssp=ss*sp;
    crot=cos(rot);
    srot=sin(rot);

    // for the x, y and z coordinates, we need
    // the coordinates of the centered, non-translated molecule

    AttractRigidbody * pLigCentered = & m_centeredligand[molIndex] ; // pointer to the centered ligand
    AttractRigidbody * pLigMoved  = & m_movedligand[molIndex] ; // pointer to the rotated/translated ligand (for forces)


    assert(shift+2 < delta.size());
    for (uint i=0; i< pLigCentered->m_activeAtoms.size(); i++)
    {
        uint atomIndex = pLigCentered->m_activeAtoms[i];

        Coord3D coords = pLigCentered->GetCoords(atomIndex);
        X = coords.x;
        Y = coords.y;
        Z = coords.z;


        xar=X*crot+Y*srot;
        yar=-X*srot+Y*crot;


        pm[0][0]=-xar*cssp-yar*cp-Z*sssp ;
        pm[1][0]=xar*cscp-yar*sp+Z*sscp ;
        pm[2][0]=0.0 ;

        pm[0][1]=-xar*sscp+Z*cscp ;
        pm[1][1]=-xar*sssp+Z*cssp ;
        pm[2][1]=-xar*cs-Z*ss ;

        pm[0][2]=yar*cscp+xar*sp ;
        pm[1][2]=yar*cssp-xar*cp ;
        pm[2][2]=-yar*ss ;


        for (uint j=0;j<3;j++)
        {
            delta[j+shift] += pm[0][j] * pLigMoved->m_forces[atomIndex].x ;
            delta[j+shift] += pm[1][j] * pLigMoved->m_forces[atomIndex].y ;
            delta[j+shift] += pm[2][j] * pLigMoved->m_forces[atomIndex].z ;
        }
    }


    if (print) std::cout << "Rotational forces: " << delta[shift] << " " << delta[shift+1] << " " << delta[shift+2] << std::endl;

    return;
}



void BaseAttractForceField::AddLigand(AttractRigidbody & lig)
{

    setDummyTypeList(lig); // sets the dummy atom type. (virtual function customized for each Attract forcefield)

    AttractRigidbody centeredlig = lig ;
    Coord3D com = lig.FindCenter();
    m_ligcenter.push_back(com);

    m_movedligand.push_back(lig);
    centeredlig.CenterToOrigin();
    m_centeredligand.push_back(centeredlig);

}



void BaseAttractForceField::MakePairLists()
{
//at this point we expect that m_movedligand still contains original coordinates of all ligands
//(ie not centered) because we will generate the pairlist from this vector (list)


//creates the pairlist: loop over all pairs of ligands
    for (uint i=0; i < m_movedligand.size(); i++)
        for (uint j=i+1; j<m_movedligand.size(); j++)
        {
            AttractPairList plist(m_movedligand[i], m_movedligand[j], m_cutoff);
            m_pairlists.push_back(plist);
        }

}


AttractRigidbody BaseAttractForceField::GetLigand(uint i) {return m_movedligand[i];};


void AttractForceField2::setDummyTypeList(AttractRigidbody& lig)
{
    lig.setDummyTypes(m_params->_dummytypes);
}


} //namespace PTools
