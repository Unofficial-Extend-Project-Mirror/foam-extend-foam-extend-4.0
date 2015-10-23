/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

InClass
    vtkPV3Foam

\*---------------------------------------------------------------------------*/

#ifndef vtkFOAMTupleRemap_H
#define vtkFOAMTupleRemap_H

// FOAM includes
#include "StaticAssert.H"
#include "Swap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline void vtkFOAMTupleRemap(float vec[]);


// a symmTensor specialization to remap FOAM -> ParaView naming order
// Qt/Core/pqScalarBarRepresentation.cxx defines this order
//     { "XX", "YY", "ZZ", "XY", "YZ", "XZ" }
// in pqScalarBarRepresentation::getDefaultComponentLabel()
// whereas FOAM uses this order
//     { XX, XY, XZ, YY, YZ, ZZ }
//
// for extra safety, assert that symmTensor indeed has 6 components
StaticAssert(Foam::symmTensor::nComponents == 6);


// Template specialization for symmTensor
template<>
inline void vtkFOAMTupleRemap<Foam::symmTensor>(float vec[])
{
    Foam::Swap(vec[1], vec[3]);   // swap XY <-> YY
    Foam::Swap(vec[2], vec[5]);   // swap XZ <-> ZZ
}


template<class Type>
inline void vtkFOAMTupleRemap(float vec[])
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
