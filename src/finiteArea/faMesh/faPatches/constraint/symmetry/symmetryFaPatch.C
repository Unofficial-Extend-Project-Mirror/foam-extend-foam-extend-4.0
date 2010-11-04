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

#include "symmetryFaPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(symmetryFaPatch, 0);
addToRunTimeSelectionTable(faPatch, symmetryFaPatch, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void symmetryFaPatch::makeCorrVecs(vectorField& cv) const
{
    // Non-orthogonal correction not allowed.  HJ, 16/Apr/2009
    cv = vector::zero;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

// Construct from components
symmetryFaPatch::symmetryFaPatch
(
    const word& name,
    const labelList& edgeLabels,
    const label index,
    const faBoundaryMesh& bm,
    const label ngbPolyPatchIndex
)
:
    faPatch(name, edgeLabels, index, bm, ngbPolyPatchIndex)
{}

//- Construct from dictionary
symmetryFaPatch::symmetryFaPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faBoundaryMesh& bm
)
:
    faPatch(name, dict, index, bm)
{
    if(ngbPolyPatchIndex() == -1)
    {
        FatalErrorIn
        (
            "symmetryFaPatch::symmetryFaPatch(const word&, const dictionary&, const label, const faBoundaryMesh&)"
        )   << "Neighbour polyPatch index is not specified for faPatch "
            << this->name() << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
