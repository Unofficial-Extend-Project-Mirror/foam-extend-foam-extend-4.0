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

Global
    Foam::preservePatchTypes

Description
    preservePatchTypes

SourceFiles
    preservePatchTypes.C

\*---------------------------------------------------------------------------*/

#ifndef preservePatchTypes_H
#define preservePatchTypes_H

#include "fileName.H"
#include "wordList.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class objectRegistry;

//- preserve patch types
void preservePatchTypes
(
    const objectRegistry& obr,
    const word& meshInstance,
    const fileName& meshDir,
    const wordList& patchNames,
    PtrList<dictionary>& patchDicts,
    const word& defaultFacesName,
    word& defaultFacesType
);

//- preserve patch types
void preservePatchTypes
(
    const objectRegistry& obr,
    const word& meshInstance,
    const fileName& meshDir,
    const wordList& patchNames,
    wordList& patchTypes,
    const word& defaultFacesName,
    word& defaultFacesType,
    wordList& patchPhysicalTypes
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
