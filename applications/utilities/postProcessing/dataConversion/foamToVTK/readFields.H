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
    Foam::readFields

Description

SourceFiles
    readFields.C

\*---------------------------------------------------------------------------*/

#ifndef readFields_H
#define readFields_H

#include "fvMesh.H"
#include "PtrList.H"
#include "IOobjectList.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Read the fields and optionally subset and put on the pointer list
template<class GeoField>
void readFields
(
    const vtkMesh& vMesh,
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields,
    PtrList<GeoField>& fields
);

// Read the fields and optionally subset and put on the pointer list
template<class GeoField>
void readFieldsNoSubset
(
    const vtkMesh& vMesh,
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields,
    PtrList<GeoField>& fields
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "readFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
