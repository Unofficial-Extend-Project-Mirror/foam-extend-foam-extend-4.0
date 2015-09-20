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

Description
    FV edge mapper.

\*---------------------------------------------------------------------------*/

#include "faEdgeMapper.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faEdgeMapper::calcAddressing() const
{
    if (directAddrPtr_)
    {
        FatalErrorIn("void faEdgeMapper::calcAddressing() const)")
            << "Addressing already calculated"
            << abort(FatalError);
    }

    // Dummy mapping: take value from edge 0
    directAddrPtr_ = new labelList(size(), 0);
}


void Foam::faEdgeMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::faEdgeMapper::faEdgeMapper
(
    const faMesh& mesh,
    const mapPolyMesh& mpm
)
:
    mesh_(mesh),
    mpm_(mpm),
    sizeBeforeMapping_(mesh.nInternalEdges()),
    directAddrPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faEdgeMapper::~faEdgeMapper()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::unallocLabelList& Foam::faEdgeMapper::directAddressing() const
{
    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}


const Foam::labelListList& Foam::faEdgeMapper::addressing() const
{
    FatalErrorIn
    (
        "const labelListList& faEdgeMapper::addressing() const"
    )   << "Requested interpolative addressing for a direct mapper."
        << abort(FatalError);

    return labelListList::null();
}


const Foam::scalarListList& Foam::faEdgeMapper::weights() const
{
    FatalErrorIn
    (
        "const scalarListList& faEdgeMapper::weights() const"
    )   << "Requested interpolative weights for a direct mapper."
        << abort(FatalError);

    return scalarListList::null();
}


// ************************************************************************* //
