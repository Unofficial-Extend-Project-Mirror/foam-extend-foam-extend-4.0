/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "passiveProcessorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "SubField.H"
#include "matchPoints.H"
#include "OFstream.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "foamTime.H"
#include "transformList.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(passiveProcessorPolyPatch, 0);

    addToRunTimeSelectionTable
    (
        polyPatch,
        passiveProcessorPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::passiveProcessorPolyPatch::passiveProcessorPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const int myProcNo,
    const int neighbProcNo,
    const labelList& globalFaceIndex
)
:
    polyPatch(name, size, start, index, bm),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo),
    globalFaceIndex_(globalFaceIndex)
{
    if (globalFaceIndex.size() != size)
    {
        FatalErrorIn
        (
            "passiveProcessorPolyPatch::passiveProcessorPolyPatch(...)"
        )   << "Bad global index list.  Patch size: " << this->size()
            << " global index size: " << globalFaceIndex.size()
            << abort(FatalError);
    }
}


Foam::passiveProcessorPolyPatch::passiveProcessorPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, dict, index, bm),
    myProcNo_(readLabel(dict.lookup("myProcNo"))),
    neighbProcNo_(readLabel(dict.lookup("neighbProcNo"))),
    globalFaceIndex_("globalFaceIndex", dict, size())
{}


Foam::passiveProcessorPolyPatch::passiveProcessorPolyPatch
(
    const passiveProcessorPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    globalFaceIndex_(newSize, -1)  // Cannot set global index.  HJ, 4/May/2018
{}


Foam::passiveProcessorPolyPatch::passiveProcessorPolyPatch
(
    const passiveProcessorPolyPatch& pp
)
:
    polyPatch(pp),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    globalFaceIndex_(pp.size(), -1)
{}


Foam::passiveProcessorPolyPatch::passiveProcessorPolyPatch
(
    const passiveProcessorPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    globalFaceIndex_(pp.size(), -1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::passiveProcessorPolyPatch::~passiveProcessorPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::passiveProcessorPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    os.writeKeyword("myProcNo") << myProcNo_
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbProcNo") << neighbProcNo_
        << token::END_STATEMENT << nl;
    globalFaceIndex_.writeEntry("globalFaceIndex", os);
}


// ************************************************************************* //
