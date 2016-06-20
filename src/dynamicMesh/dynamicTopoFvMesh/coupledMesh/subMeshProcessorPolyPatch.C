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

Class
    subMeshProcessorPolyPatch

Description
    Functions specific to subMeshProcessorPolyPatch

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "subMeshProcessorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(subMeshProcessorPolyPatch, 0);
addToRunTimeSelectionTable(polyPatch, subMeshProcessorPolyPatch, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

subMeshProcessorPolyPatch::subMeshProcessorPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const int myProcNo,
    const int neighbProcNo
)
:
    processorPolyPatch
    (
        name,
        size,
        start,
        index,
        bm,
        myProcNo,
        neighbProcNo
    )
{}


subMeshProcessorPolyPatch::subMeshProcessorPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    processorPolyPatch(name, dict, index, bm)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

subMeshProcessorPolyPatch::~subMeshProcessorPolyPatch()
{}

} // End namespace Foam

// ************************************************************************* //
