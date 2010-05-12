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

Description

\*---------------------------------------------------------------------------*/

#include "UILList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
UILList<LListBase, T>::UILList(const UILList<LListBase, T>& slpl)
{
    for (const_iterator iter = slpl.begin(); iter != slpl.end(); ++iter)
    {
        append(&iter());
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class LListBase, class T>
void UILList<LListBase, T>::operator=(const UILList<LListBase, T>& slpl)
{
    LListBase::clear();

    for (const_iterator iter = slpl.begin(); iter != slpl.end(); ++iter)
    {
        append(&iter());
    }
}


// Comparison for equality
template<class LListBase, class T>
bool UILList<LListBase, T>::operator==(const UILList<LListBase, T>& slpl) const
{
    if (this->size() != slpl.size())
    {
        return false;
    }

    bool equal = true;

    const_iterator iter1 = this->begin();
    const_iterator iter2 = slpl.begin();

    for (; iter1 != this->end(); ++iter1, ++iter2)
    {
        equal = equal && iter1() == iter2();
    }

    return equal;
}


// Comparison for inequality
template<class LListBase, class T>
bool UILList<LListBase, T>::operator!=(const UILList<LListBase, T>& slpl) const
{
    return !operator==(slpl);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "UILListIO.C"


// ************************************************************************* //
