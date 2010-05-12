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
    PolyPatch mapper for the face tetFem decomposition

\*---------------------------------------------------------------------------*/

#include "tetPolyPatchMapperFaceDecomp.H"
#include "tetPolyPatchFaceDecomp.H"
#include "tetPolyBoundaryMeshFaceDecomp.H"
#include "tetPolyMeshFaceDecomp.H"
#include "mapPolyMesh.H"
#include "pointMapper.H"
#include "faceMapper.H"
#include "faceTetPolyPatchFaceDecomp.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::tetPolyPatchMapperFaceDecomp::calcAddressing() const
{
    if
    (
        directPtr_
     || directAddrPtr_
     || interpolationAddrPtr_
     || weightsPtr_
    )
    {
        FatalErrorIn
        (
            "void tetPolyPatchMapperFaceDecomp::calcAddressing() const)"
        )   << "Addressing already calculated"
            << abort(FatalError);
    }

    // Mapping

    // Calculate direct (if all are direct)
    directPtr_ = new bool(pMapper_.direct() && fMapper_.direct());

    const labelList& curPatchPointMap = mpm_.patchPointMap()[patch_.index()];

    const label patchOffset = mpm_.oldPatchStarts()[patch_.index()];
    const label oldPatchFaceOffset = mpm_.oldPatchNMeshPoints()[patch_.index()];

    // Assemble the maps
    // If it's a face patch, insert the faces
    const polyPatch& p =
        refCast<const faceTetPolyPatchFaceDecomp>(patch_).patch();

    const label curPatchStart = p.start();
    const label curPatchEnd = curPatchStart + p.size();

    // Mapping
    const label oldPatchStart = mpm_.oldPatchStarts()[patch_.index()];

    const label oldPatchEnd =
        oldPatchStart + mpm_.oldPatchSizes()[patch_.index()];

    if (*directPtr_)
    {
        // Direct mapping

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

        const labelList& mappedFaces = fMapper_.directAddressing();

        // Insert faces
        for (label faceI = curPatchStart; faceI < curPatchEnd; faceI++)
        {
            if
            (
                mappedFaces[faceI] >= oldPatchStart
             && mappedFaces[faceI] < oldPatchEnd
            )
            {
                addr[nAddr] =
                    mappedFaces[faceI] - patchOffset + oldPatchFaceOffset;
            }
            else
            {
                addr[nAddr] = 0;
            }

            nAddr++;
        }
    }
    else
    {
        // Interpolative mapping
        interpolationAddrPtr_ = new labelListList(size());
        labelListList& addr = *interpolationAddrPtr_;

        weightsPtr_ = new scalarListList(size());
        scalarListList& w = *weightsPtr_;

        label nAddr = 0;

        // Insert points

        forAll (curPatchPointMap, pointI)
        {
            if (curPatchPointMap[pointI] > -1)
            {
                addr[nAddr] = labelList(1, curPatchPointMap[pointI]);
            }
            else
            {
                addr[nAddr] = labelList(1, 0);
            }

            w[nAddr] = scalarList(1, 1.0);
            nAddr++;
        }

        const labelListList& mappedFaces = fMapper_.addressing();
        const scalarListList& faceWeights = fMapper_.weights();

        // Insert faces
        for (label faceI = curPatchStart; faceI < curPatchEnd; faceI++)
        {
            labelList& curAddr = addr[nAddr];
            scalarList& curW = w[nAddr];
            label nActive = 0;

            const labelList& curMf = mappedFaces[faceI];
            curAddr.setSize(curMf.size());

            const scalarList& curMfWeights = faceWeights[faceI];
            curW.setSize(curMfWeights.size());

            forAll (curMf, lfI)
            {
                if
                (
                    curMf[lfI] >= oldPatchStart
                 && curMf[lfI] < oldPatchEnd
                )
                {
                    curAddr[nActive] =
                        curMf[lfI] - patchOffset + oldPatchFaceOffset;

                    curW[nActive] = curMfWeights[lfI];
                    nActive++;
                }

            }

            // Cater for bad mapping
            if (nActive == 0)
            {
                curAddr[nActive] = 0;
                curW[nActive] = 1;
                nActive++;
            }

            curAddr.setSize(nActive);
            curW.setSize(nActive);

            // Re-scale the weights
            scalar sumW = sum(curW);

            forAll (curW, wI)
            {
                curW[wI] /= sumW;
            }

            nAddr++;
        }
    }
}


void Foam::tetPolyPatchMapperFaceDecomp::clearOut()
{
    deleteDemandDrivenData(directPtr_);
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::tetPolyPatchMapperFaceDecomp::tetPolyPatchMapperFaceDecomp
(
    const tetPolyPatchFaceDecomp& patch,
    const mapPolyMesh& meshMap,
    const pointMapper& pMapper,
    const faceMapper& fMapper
)
:
    patch_(patch),
    mpm_(meshMap),
    pMapper_(pMapper),
    fMapper_(fMapper),
    directPtr_(NULL),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tetPolyPatchMapperFaceDecomp::~tetPolyPatchMapperFaceDecomp()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::tetPolyPatchMapperFaceDecomp::size() const
{
    return patch_.size();
}


Foam::label Foam::tetPolyPatchMapperFaceDecomp::sizeBeforeMapping() const
{
    return
        mpm_.oldPatchSizes()[patch_.index()]
      + mpm_.oldPatchNMeshPoints()[patch_.index()];
}


bool Foam::tetPolyPatchMapperFaceDecomp::direct() const
{
    if (!directPtr_)
    {
        calcAddressing();
    }

    return *directPtr_;
}


const Foam::unallocLabelList&
Foam::tetPolyPatchMapperFaceDecomp::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const unallocLabelList& tetPolyPatchMapperFaceDecomp::"
            "directAddressing() const"
        )   << "Requested direct addressing for an interpolative mapper."
            << abort(FatalError);
    }

    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}


const Foam::labelListList&
Foam::tetPolyPatchMapperFaceDecomp::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& tetPolyPatchMapperFaceDecomp::"
            "addressing() const"
        )   << "Requested interpolative addressing for a direct mapper."
            << abort(FatalError);
    }

    if (!interpolationAddrPtr_)
    {
        calcAddressing();
    }

    return *interpolationAddrPtr_;
}


const Foam::scalarListList& Foam::tetPolyPatchMapperFaceDecomp::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& tetPolyPatchMapperFaceDecomp::"
            "weights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcAddressing();
    }

    return *weightsPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
