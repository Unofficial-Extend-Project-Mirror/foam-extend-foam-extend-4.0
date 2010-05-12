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
    PolyPatch mapper for the cell tetFem decomposition

\*---------------------------------------------------------------------------*/

#include "tetPolyPatchMapperCellDecomp.H"
#include "tetPolyPatchCellDecomp.H"
#include "tetPolyBoundaryMeshCellDecomp.H"
#include "tetPolyMeshCellDecomp.H"
#include "tetFemMatrices.H"
#include "mapPolyMesh.H"
#include "pointMapper.H"
#include "faceTetPolyPatchCellDecomp.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::tetPolyPatchMapperCellDecomp::calcAddressing() const
{
    if (directAddrPtr_)
    {
        FatalErrorIn
        (
            "void tetPolyPatchMapperCellDecomp::calcAddressing() const)"
        )   << "Addressing already calculated"
            << abort(FatalError);
    }

    // Mapping

    const labelList& curPatchPointMap = mpm_.patchPointMap()[patch_.index()];

    // Assemble the map (direct mapping)
    directAddrPtr_ = new labelList(size());
    labelList& addr = *directAddrPtr_;
    label nAddr = 0;

    forAll (curPatchPointMap, pointI)
    {
        if (curPatchPointMap[pointI] > -1)
        {
            addr[nAddr] = curPatchPointMap[pointI];
        }
        else
        {
            addr[nAddr] = 0;
        }
        nAddr++;
    }
}


void Foam::tetPolyPatchMapperCellDecomp::clearOut()
{
    deleteDemandDrivenData(directPtr_);
    deleteDemandDrivenData(directAddrPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::tetPolyPatchMapperCellDecomp::tetPolyPatchMapperCellDecomp
(
    const tetPolyPatchCellDecomp& patch,
    const mapPolyMesh& meshMap
)
:
    patch_(patch),
    mpm_(meshMap),
    directPtr_(NULL),
    directAddrPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tetPolyPatchMapperCellDecomp::~tetPolyPatchMapperCellDecomp()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::tetPolyPatchMapperCellDecomp::size() const
{
    return patch_.size();
}


Foam::label Foam::tetPolyPatchMapperCellDecomp::sizeBeforeMapping() const
{
    return mpm_.oldPatchNMeshPoints()[patch_.index()];
}


const Foam::unallocLabelList&
Foam::tetPolyPatchMapperCellDecomp::directAddressing() const
{
    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
