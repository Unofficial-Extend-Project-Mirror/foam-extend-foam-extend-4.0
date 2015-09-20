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

InApplication
    foamToEnsight

Description

SourceFiles
    ensightField.C

\*---------------------------------------------------------------------------*/

#ifndef ensightField_H
#define ensightField_H

#include "ensightMesh.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void ensightField
(
    const Foam::IOobject& fieldObject,
    const Foam::ensightMesh& eMesh,
    const Foam::fileName& postProcPath,
    const Foam::word& prepend,
    const Foam::label timeIndex,
    const bool binary,
    Foam::Ostream& ensightCaseFile
);


template<class Type>
void writePatchField
(
    const Foam::word& fieldName,
    const Foam::Field<Type>& pf,
    const Foam::word& patchName,
    const Foam::ensightMesh& eMesh,
    const Foam::fileName& postProcPath,
    const Foam::word& prepend,
    const Foam::label timeIndex,
    Foam::Ostream& ensightCaseFile
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ensightField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
