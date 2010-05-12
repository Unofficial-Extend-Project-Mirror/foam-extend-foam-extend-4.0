/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "LPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
LPtrList<LListBase, T>::LPtrList(const LPtrList<LListBase, T>& slpl)
:
    LList<LListBase, T*>()
{
    for(const_iterator iter = slpl.begin(); iter != slpl.end(); ++iter)
    {
        append(iter().clone().ptr());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class LListBase, class T>
LPtrList<LListBase, T>::~LPtrList()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return and remove head
template<class LListBase, class T>
bool LPtrList<LListBase, T>::eraseHead()
{
    T* tPtr;
    if ((tPtr = this->removeHead()))
    {
        delete tPtr;
        return true;
    }
    else
    {
        return false;
    }
}


template<class LListBase, class T>
void LPtrList<LListBase, T>::clear()
{
    label oldSize = this->size();
    for (label i=0; i<oldSize; i++)
    {
        eraseHead();
    }

    LList<LListBase, T*>::clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class LListBase, class T>
void LPtrList<LListBase, T>::operator=(const LPtrList<LListBase, T>& slpl)
{
    clear();

    for(const_iterator iter = slpl.begin(); iter != slpl.end(); ++iter)
    {
        append(iter().clone().ptr());
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "LPtrListIO.C"


// ************************************************************************* //
