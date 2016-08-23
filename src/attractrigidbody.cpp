/****************************************************************************
 *   Copyright (C) 2008   Adrien Saladin                                    *
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

#include <iostream>
#include <sstream>

#include "attractrigidbody.h"

namespace PTools{

///////////// AttractRigidbody implementation:

AttractRigidbody::AttractRigidbody(const Rigidbody & rig)
        : Rigidbody(rig)
{
this->init_();
}


AttractRigidbody::AttractRigidbody(const std::string& filename)
 : Rigidbody(filename)
{
this->init_();
}

void AttractRigidbody::init_()
{
    // extracts the "extra" field of Atoms to the m_atomTypeNumber array:
    uint   atcategory  = 0;
    dbl  atcharge   = 0.0;

    for (uint i = 0; i < Size() ; ++i)
    {
        Atomproperty & at (mAtomProp[i]);
        std::string extra = at.extra;
        std::istringstream iss( extra );
        iss >> atcategory >> atcharge ;
        m_atomTypeNumber.push_back(atcategory-1);  // -1 to directly use the atomTypeNumber into C-array
        m_charge.push_back(atcharge);

    }

    setRotation(true);
    setTranslation(true);

    resetForces();
}


void AttractRigidbody::setDummyTypes(const std::vector<uint>& dummy)
{
    m_dummytypes = dummy;
    this->updateActiveList();
}



void AttractRigidbody::updateActiveList()
{
 std::vector<uint> newactivelist;

 for(uint i=0; i<this->Size(); i++)
 {
   if( isAtomActive(i) ) newactivelist.push_back(i);
 }

std::swap(m_activeAtoms, newactivelist);

}

void AttractRigidbody::AddAtom(const Atom& at)
{
    Atomproperty atp(at);
    Coord3D co = at.coords;
    AddAtom(atp,co);
}

void AttractRigidbody::AddAtom(const Atomproperty& at, Coord3D co)
{
    mAtomProp.push_back(at);
    AddCoord(co);
    //Add AttractRigidbody parameters
    uint   atcategory  = 0;
    dbl  atcharge   = 0.0;
    std::string extra = at.extra;
    std::istringstream iss( extra );
    iss >> atcategory >> atcharge ;
    m_atomTypeNumber.push_back(atcategory-1);  // -1 to directly use the atomTypeNumber into C-array
    m_charge.push_back(atcharge);
    m_forces.push_back(Coord3D());
}

} //namespace PTools


