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
    Point mapper for the face tetFem decomposition

\*---------------------------------------------------------------------------*/

#include "tetPointMapperFaceDecomp.H"
#include "tetPolyMeshFaceDecomp.H"
#include "mapPolyMesh.H"
#include "pointMapper.H"
#include "faceMapper.H"
#include "cellMapper.H"
#include "tetPointMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::tetPointMapperFaceDecomp::calcAddressing() const
{
    if
    (
        directPtr_
     || directAddrPtr_
     || interpolationAddrPtr_
     || weightsPtr_
     || insertedObjectsPtr_
     || insertedObjectLabelsPtr_
    )
    {
        FatalErrorIn("void tetPointMapperFaceDecomp::calcAddressing() const)")
            << "Addressing already calculated"
            << abort(FatalError);
    }

    const label oldFaceOffset = mpm_.nOldPoints();
    const label oldCellOffset = oldFaceOffset + mpm_.nOldFaces();

    // Mapping

    // Calculate direct (if all are direct)
    directPtr_ =
        new bool
        (
            pointMap_.direct()
         && faceMap_.direct()
         && cellMap_.direct()
        );

    // Assemble the maps
    if (*directPtr_)
    {
        // Direct mapping
        const labelList& mappedPoints = pointMap_.directAddressing();
        const labelList& mappedFaces = faceMap_.directAddressing();
        const labelList& mappedCells = cellMap_.directAddressing();

        directAddrPtr_ = new labelList(size());
        labelList& addr = *directAddrPtr_;
        label nAddr = 0;

        forAll (mappedPoints, pointI)
        {
            addr[nAddr] = mappedPoints[pointI];
            nAddr++;
        }

        forAll (mappedFaces, faceI)
        {
            addr[nAddr] = mappedFaces[faceI] + oldFaceOffset;
            nAddr++;
        }

        forAll (mappedCells, cellI)
        {
            addr[nAddr] = mappedCells[cellI] + oldCellOffset;
            nAddr++;
        }
    }
    else
    {
        // Interpolative mapping
        if (tetPolyMeshFaceDecomp::debug)
        {
            Info<< " interpolative" << endl;
        }

        // Interpolative mapping
        interpolationAddrPtr_ = new labelListList(size());
        labelListList& addr = *interpolationAddrPtr_;

        weightsPtr_ = new scalarListList(size());
        scalarListList& w = *weightsPtr_;

        label nAdded = 0;

        // Insert points
        const labelList& mappedPoints = pointMap_.directAddressing();

        forAll (mappedPoints, pointI)
        {
            addr[nAdded] = labelList(1, mappedPoints[pointI]);
            w[nAdded] = scalarList(1, 1.0);
            nAdded++;
        }

        // Do face addressing, direct or interpolative
        if (faceMap_.direct())
        {
            // Direct for faces
            const labelList& mappedFaces = faceMap_.directAddressing();

            // Insert faces
            forAll (mappedFaces, faceI)
            {
                addr[nAdded] =
                    labelList(1, mappedFaces[faceI]  + oldFaceOffset);
                w[nAdded] = scalarList(1, 1.0);
                nAdded++;
            }
        }
        else
        {

            // Interpolative for faces
            const labelListList& mappedFaces = faceMap_.addressing();
            const scalarListList& faceWeights = faceMap_.weights();

            // Insert faces
            forAll (mappedFaces, faceI)
            {
                labelList& curAddr = addr[nAdded];

                const labelList& curMf = mappedFaces[faceI];
                curAddr.setSize(curMf.size());

                forAll (curAddr, cI)
                {
                    curAddr[cI] = curMf[cI] + oldFaceOffset;
                }

                // Weights remain the same
                w[nAdded] = faceWeights[faceI];
                nAdded++;
            }
        }

        // Do cell addressing, direct or interpolative
        if (cellMap_.direct())
        {
            // Direct for cells
            const labelList& mappedCells = cellMap_.directAddressing();

            // Insert cells
            forAll (mappedCells, cellI)
            {
                addr[nAdded] =
                    labelList(1, mappedCells[cellI] + oldCellOffset);
                w[nAdded] = scalarList(1, 1.0);
                nAdded++;
            }
        }
        else
        {
            // Interpolative for cells
            const labelListList& mappedCells = cellMap_.addressing();
            const scalarListList& cellWeights = cellMap_.weights();

            // Insert cells
            forAll (mappedCells, cellI)
            {
                labelList& curAddr = addr[nAdded];

                const labelList& curMc = mappedCells[cellI];
                curAddr.setSize(curMc.size());

                forAll (curAddr, cI)
                {
                    curAddr[cI] = curMc[cI] + oldCellOffset;
                }

                // Weights remain the same
                w[nAdded] = cellWeights[cellI];
                nAdded++;
            }
        }
    }

    // Inserted objects

    // Are there inserted objects presents
    insertedObjectsPtr_ =
        new bool
        (
            pointMap_.insertedObjects()
         || faceMap_.insertedObjects()
         || cellMap_.insertedObjects()
        );

    // If there are, assemble the labels
    if (*insertedObjectsPtr_)
    {
        const labelList& insPoints = pointMap_.insertedObjectLabels();
        const labelList& insFaces = faceMap_.insertedObjectLabels();
        const labelList& insCells = cellMap_.insertedObjectLabels();

        insertedObjectLabelsPtr_ =
            new labelList(insPoints.size() + insFaces.size() + insCells.size());
        labelList& ins = *insertedObjectLabelsPtr_;

        label nIns = 0;

        forAll (insPoints, pointI)
        {
            ins[nIns] = insPoints[pointI];
            nIns++;
        }

        forAll (insFaces, faceI)
        {
            ins[nIns] = insFaces[faceI] + oldFaceOffset;
            nIns++;
        }

        forAll (insCells, cellI)
        {
            ins[nIns] = insCells[cellI] + oldCellOffset;
            nIns++;
        }
    }
    else
    {
        // No inserted objects
        insertedObjectLabelsPtr_ = new labelList(0);
    }
}


void Foam::tetPointMapperFaceDecomp::clearOut()
{
    deleteDemandDrivenData(directPtr_);
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);

    deleteDemandDrivenData(insertedObjectsPtr_);
    deleteDemandDrivenData(insertedObjectLabelsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::tetPointMapperFaceDecomp::tetPointMapperFaceDecomp
(
    const tetPolyMeshFaceDecomp& mesh,
    const mapPolyMesh& meshMap,
    const pointMapper& pMapper,
    const faceMapper& fMapper,
    const cellMapper& cMapper
)
:
    mesh_(mesh),
    mpm_(meshMap),
    pointMap_(pMapper),
    faceMap_(fMapper),
    cellMap_(cMapper),
    size_(mesh().nPoints() + mesh().nFaces() + mesh().nCells()),
    directPtr_(NULL),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    insertedObjectsPtr_(NULL),
    insertedObjectLabelsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tetPointMapperFaceDecomp::~tetPointMapperFaceDecomp()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::tetPointMapperFaceDecomp::size() const
{
    return size_;
}


Foam::label Foam::tetPointMapperFaceDecomp::sizeBeforeMapping() const
{
    return mpm_.nOldPoints() + mpm_.nOldFaces() + mpm_.nOldCells();
}


bool Foam::tetPointMapperFaceDecomp::direct() const
{
    if (!directPtr_)
    {
        calcAddressing();
    }

    return *directPtr_;
}


const Foam::unallocLabelList&
Foam::tetPointMapperFaceDecomp::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const unallocLabelList& tetPointMapperFaceDecomp::"
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


const Foam::labelListList& Foam::tetPointMapperFaceDecomp::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& tetPointMapperFaceDecomp::addressing() const"
        )   << "Requested interpolative addressing for a direct mapper."
            << abort(FatalError);
    }

    if (!interpolationAddrPtr_)
    {
        calcAddressing();
    }

    return *interpolationAddrPtr_;
}


const Foam::scalarListList& Foam::tetPointMapperFaceDecomp::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& tetPointMapperFaceDecomp::weights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcAddressing();
    }

    return *weightsPtr_;
}


bool Foam::tetPointMapperFaceDecomp::insertedObjects() const
{
    if (!insertedObjectsPtr_)
    {
        calcAddressing();
    }

    return *insertedObjectsPtr_;
}


const Foam::labelList&
Foam::tetPointMapperFaceDecomp::insertedObjectLabels() const
{
    if (!insertedObjectLabelsPtr_)
    {
        calcAddressing();
    }

    return *insertedObjectLabelsPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
