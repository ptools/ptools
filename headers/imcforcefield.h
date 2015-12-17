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


#ifndef IMCFORCEFIELD_H
#define IMCFORCEFIELD_H

#include "attractforcefield.h"


namespace PTools {



class ImcForceField: public BaseAttractForceField
{



public:

    ImcForceField(std::string paramsFileName, dbl cutoff); //constructor
    void InitParams(const std::string & paramsFileName); //read and initalize forcefield parameters

    dbl nonbon8_forces(AttractRigidbody& rec, AttractRigidbody& lig, AttractPairList & pairlist, std::vector<Coord3D>& forcerec, std::vector<Coord3D>& forcelig, bool print=false);





    virtual ~ImcForceField() {};

private:

    struct Parameters // all parameters needed for bead-bead interaction
    {
        uint id1;  //id of first atom
        uint id2;
        int Q;
        dbl A;
        uint N;   // number of gaussians (redundant with B.size())

        std::vector<dbl> B;
        std::vector<dbl> X;
        std::vector<dbl> C;

    };

    dbl m_rstk;
    
    void setDummyTypeList(AttractRigidbody& lig) {
        std::vector<uint> dummytypes;
        lig.setDummyTypes(dummytypes);
    }; //forcefield1 has no dummy type

    Parameters m_parameters[33][33] ;
    int m_nbparameters[33][33]; // number of gaussians needed to describe the interaction between two atom types



};



} // namespace PTools



#endif //IMCFORCEFIELD_H
