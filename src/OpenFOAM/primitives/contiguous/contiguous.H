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

InClass
    Foam::contiguous

Description
    Template function which specifies if the data of the type is contiguous
    or not.

    The default function specifies that data are not contiguous.
    This is specialised for the types (e.g. primitives) for which their
    data are contiguous.

\*---------------------------------------------------------------------------*/

#ifndef contiguous_H
#define contiguous_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Assume the data associated with type T is not contiguous
template<class T>
inline bool contiguous()              {return false;}


// Specify data associated with primitive types is contiguous

template<>
inline bool contiguous<bool>()        {return true;}

template<>
inline bool contiguous<char>()        {return true;}

template<>
inline bool contiguous<short>()       {return true;}

template<>
inline bool contiguous<int>()         {return true;}

template<>
inline bool contiguous<long>()        {return true;}

template<>
inline bool contiguous<float>()       {return true;}

template<>
inline bool contiguous<double>()      {return true;}

template<>
inline bool contiguous<long double>() {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
