/****************************************************************************
 *   Copyright (C) 2014   Adrien Saladin                                    *
 *   adrien.saladin@gmail.com                                               *
 *                                                                          *
 *   This program is free software: you can redistribute it and/or modify   *
 *   it under the terms of the GNU General Public License as published by   *
 *   the Free Software Foundation, either version 3 of the License, or      *
 *   (at your option) any later version.                                    *
 *                                                                          *
 *   This program is distributed in the hope that it will be useful,        *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU General Public License for more details.                           *
 *                                                                          *
 *   You should have received a copy of the GNU General Public License      *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                          *
 ***************************************************************************/
/*
 *
 * Implementation of IMC/US force field
 *
 *  Naome A. 2013
 *
 *  Functional form:
 *  U(x) = A/x**12 + B1*exp(-C1*(x-x1)**2) + B2*exp(-C2*(x-x2)**2) + B3*exp(-C3*(x-x3)**2) + B4*exp(-C4*(x-x4)**2) + B5*exp(-C5*(x-x5)**2) + Q*138.935485/(80*x)
 */

#include "imcforcefield.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>



using namespace std;

namespace PTools
{





ImcForceField::ImcForceField(std::string paramsFileName, dbl cutoff)
{

    InitParams(paramsFileName);
    
    m_rstk=0.0; //no restraint by default
    m_cutoff = cutoff;

}


void ImcForceField::InitParams(const std::string & paramsFileName )
{

    std::ifstream aminon(paramsFileName.c_str());
    if (!aminon)
    {
        //the file cannot be opened
        std::string msg = "imcforcefield.cpp: Cannot Locate file forcefield parameters (imc.par)\n";
        std::cerr << msg ;
        throw std::invalid_argument(msg);
    }

    uint count = 0;
    uint linenb = 0;

    while (!aminon.eof())
    {

        std::string line;
        std::getline(aminon, line);
        //cout << line << endl;

        linenb++;

        if (line[0] == '#')
        {
//            cout << "skipping line " << linenb << endl;
            continue;
        }


        count ++;

        Parameters p;
        std::string type1, type2;
        dbl Ai,Bi,xi,Ci;

        std::istringstream linestream(line);

        linestream >> p.id1 >> p.id2 >> type1 >> type2;
        linestream >> p.Q >> Ai >> p.N;
	p.A = Ai*1e12/4.184; // from kJ.mol^-1.nm^12 to kcal.mol^-1.A^12

        p.id1--;
        p.id2--; //convert ids to C-style

//         cout <<"at line " << linenb << " id1, id2: " << p.id1 <<" "<< p.id2 << endl;

        //assert(p.id1 < 32);
        //assert(p.id2 < 32);

         //cout << "**id: " << p.id1 << " "  << p.id2 << ";Q: " << p.Q << "; A="<< p.A <<  "; N=" << p.N << " "<<endl;

        while(!linestream.eof())
        {

            linestream >> Bi >> xi >> Ci; // ! parameters are in kJ and nm units !
            p.B.push_back(Bi/4.184); // from kJ.mol^-1 to kcal.mol^-1
            p.X.push_back(xi*10.0); // from nm to A
            p.C.push_back(Ci/100.0); // from nm^-2 to A^-2

          //cout << Bi << " " << xi << " " << Ci << " " <<endl;
        }

//         cout << endl;


        this->m_parameters[p.id1][p.id2] = p;
        this->m_parameters[p.id2][p.id1] = p;  //this is a symetric array
	if ((p.id1==32) && (p.id2==32)) break;

    }

}


dbl ImcForceField::nonbon8_forces(AttractRigidbody& rec, AttractRigidbody& lig, AttractPairList & pairlist, std::vector<Coord3D>& forcerec, std::vector<Coord3D>& forcelig, bool print)
{

    assert(forcerec.size() == rec.Size());
    assert(forcelig.size() == lig.Size());

    dbl sumLJ=0.0 ; // initialize total vdW energy
    dbl sumElectrostatic=0.0; // initialize total electrostatic energy


    //synchronize coordinates for using unsafeGetCoords
    rec.syncCoords();
    lig.syncCoords();

    Coord3D a, b;
    
    //cout << pairlist.Size() << " pairs"<< endl;
    for (uint iter=0; iter<pairlist.Size(); iter++) // loop over all pairs
    {

        uint ir = pairlist[iter].atrec;
        uint jl = pairlist[iter].atlig;

        uint rAtomCat = rec.getAtomTypeNumber(ir);
        uint lAtomCat = lig.getAtomTypeNumber(jl);

	uint Qij = m_parameters[rAtomCat][lAtomCat].Q; // return charge product for pairtype 
	dbl Aij = m_parameters[rAtomCat][lAtomCat].A; // return 1/r^12 parameter
	std::vector<dbl> Bij =  m_parameters[rAtomCat][lAtomCat].B; // return Gaussian parameters
	std::vector<dbl> Xij =  m_parameters[rAtomCat][lAtomCat].X;
	std::vector<dbl> Cij =  m_parameters[rAtomCat][lAtomCat].C;	
	
	//cout << "Aij: " << Aij << ", Bij: "<<  Bij[0] << endl;
		
        lig.unsafeGetCoords(jl,a);
        rec.unsafeGetCoords(ir,b);

        Coord3D dx = a-b ; // distance vector
        dbl r2 = Norm2(dx); // square of the norm of dx
    	dbl r1 = Norm(dx); // norm of dx
        
        if (r2 < 0.001 ) r2=0.001;
        //cout << "id1:"<< rAtomCat <<", id2:"<< lAtomCat << ", d:" << r1 << endl;
        // vdW part
        dbl rr1 = 1.0/r1;
        dbl rr2 = rr1*rr1; // 1/r^2       
        dbl rr26 = rr2*rr2*rr2*rr2*rr2*rr2; // 1/r^12
	dbl repwall = Aij*rr26; // repulsive wall energy
	dbl G[5]; // declare array of 5 Gaussian energies
	dbl vlj=0.0; // declare and initialize the vdW energy for current pair
	dbl fb=0.0; // declare and initialize the minus derivative of the energy
	//cout << "repwall: " << repwall << endl;
	if ((r1 <= 3.2) && (Aij!=0.0)) //at short distance some gaussians are so negative that an unphysical energy well appears, this makes sure only the repulsive wall contribution is calculated 
	{
		//cout << "id1:"<< rAtomCat <<", id2:"<< lAtomCat << ", d:" << r1 << ", repwall: " << repwall << " www " << endl;
		//continue;
	}
	else
	{
		for ( uint i=0; i<Bij.size(); i++) // loop over Gaussians
		{
			dbl gdist = r1-Xij[i];
			G[i] = Bij[i]*exp(-Cij[i]*(gdist*gdist)); // calculate energy of Gaussian i
			//cout << i << ". G: "<< G[i] << ", X: " << Xij[i] << endl;
			vlj += G[i];
			fb += 2.0*Cij[i]*Cij[i]*gdist*G[i]; // -dG/dr=2BijCij(rij-Xij)G
		}
	
        	//cout << "id1:"<< rAtomCat <<", id2:"<< lAtomCat << ", sz:" <<Bij.size()<<", d:" << r1 << ", repwall: " << repwall << ", vlj: "<<  vlj << endl;
	}
	//cout <<	"d: " << r1 << endl;

        vlj += repwall;
	sumLJ += vlj;
		
	fb += 12.0*repwall*rr1;
		
	dx = rr1*dx; // distance vector normalized to unity (dx/Norm(dx)) 
		
	Coord3D fdb = fb*dx ;

        //assign force to the atoms:
        forcelig[jl] -= fdb ;
        forcerec[ir] += fdb ;

		
	// electrostatic part
	if (fabs(Qij) > 0.0)
	{
		dbl et = Qij*(332.053986/80.0)*rr1;
		sumElectrostatic += et;
		
		Coord3D fdb = et*rr1*dx;
	               
                forcelig[jl] -= fdb ;
                forcerec[ir] += fdb ;			
	}
		
     }
     m_vdw = sumLJ;
     m_elec = sumElectrostatic;
     //cout << "energy: " << sumLJ << endl;
     return sumLJ + sumElectrostatic;

}



} //namespace PTools
