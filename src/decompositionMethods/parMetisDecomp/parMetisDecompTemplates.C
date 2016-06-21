/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "parMetisDecomp.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Insert at front of list
template<class Type>
void Foam::parMetisDecomp::prepend
(
    const UList<Type>& extraLst,
    List<Type>& lst
)
{
    label nExtra = extraLst.size();

    // Make space for initial elements
    lst.setSize(lst.size() + nExtra);
    for (label i = lst.size()-1; i >= nExtra; i--)
    {
        lst[i] = lst[i-nExtra];
    }

    // Insert at front
    forAll(extraLst, i)
    {
        lst[i] = extraLst[i];
    }
}

// Insert at back of list
template<class Type>
void Foam::parMetisDecomp::append
(
    const UList<Type>& extraLst,
    List<Type>& lst
)
{
    label sz = lst.size();

    // Make space for initial elements
    lst.setSize(sz + extraLst.size());

    // Insert at back
    forAll(extraLst, i)
    {
        lst[sz++] = extraLst[i];
    }
}


// ************************************************************************* //
