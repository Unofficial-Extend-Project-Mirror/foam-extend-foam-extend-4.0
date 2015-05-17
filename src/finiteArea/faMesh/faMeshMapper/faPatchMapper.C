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

\*---------------------------------------------------------------------------*/

#include "faPatchMapper.H"
#include "faPatch.H"
#include "faBoundaryMesh.H"
#include "faMesh.H"
#include "mapPolyMesh.H"
#include "faceMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faPatchMapper::calcAddressing() const
{
    if (directAddrPtr_)
    {
        FatalErrorIn
        (
            "void faPatchMapper::calcAddressing() const)"
        )   << "Addressing already calculated"
            << abort(FatalError);
    }

    directAddrPtr_ = new labelList(patch_.size(), 0);
    labelList& addr = *directAddrPtr_;

    // Make a map of old edgeFaces, giving edge index in patch given the new
    // face label next to the patch

    // Create edge index lookup
    Map<label> edgeIndexLookup;

    const labelList& reverseFaceMap = mpm_.reverseFaceMap();

    forAll (oldEdgeFaces_, oefI)
    {
        if (reverseFaceMap[oldEdgeFaces_[oefI]] > -1)
        {
            // Face has survived.  Insert its label under new face index
            edgeIndexLookup.insert(reverseFaceMap[oldEdgeFaces_[oefI]], oefI);
        }
    }

    // Go through new edgeFaces and for each edge try to locate old index
    const labelList& ef = patch_.edgeFaces();

    forAll (ef, efI)
    {
        if (edgeIndexLookup.found(ef[efI]))
        {
            addr[efI] = edgeIndexLookup[ef[efI]];
        }
        else
        {
            // Not found: map from zero
            addr[efI] = 0;
        }
    }
}


void Foam::faPatchMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::faPatchMapper::faPatchMapper
(
    const faPatch& patch,
    const mapPolyMesh& mpm
)
:
    patch_(patch),
    mpm_(mpm),
    sizeBeforeMapping_(patch.size()),
    oldEdgeFaces_(patch.edgeFaces()),
    directAddrPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faPatchMapper::~faPatchMapper()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::unallocLabelList& Foam::faPatchMapper::directAddressing() const
{
    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}


const Foam::labelListList& Foam::faPatchMapper::addressing() const
{
    FatalErrorIn
    (
        "const labelListList& faPatchMapper::addressing() const"
    )   << "Requested interpolative addressing for a direct mapper."
        << abort(FatalError);

    return labelListList::null();
}


const Foam::scalarListList& Foam::faPatchMapper::weights() const
{
    FatalErrorIn
    (
        "const scalarListList& faPatchMapper::weights() const"
    )   << "Requested interpolative weights for a direct mapper."
        << abort(FatalError);

    return scalarListList::null();
}


// ************************************************************************* //
