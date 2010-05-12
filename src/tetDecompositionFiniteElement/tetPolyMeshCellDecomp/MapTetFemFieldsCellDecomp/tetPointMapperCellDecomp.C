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

#include "tetPointMapperCellDecomp.H"
#include "tetPolyMeshCellDecomp.H"
#include "mapPolyMesh.H"
#include "pointMapper.H"
#include "faceMapper.H"
#include "cellMapper.H"
#include "tetPointMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::tetPointMapperCellDecomp::calcAddressing() const
{
    if (tetPolyMeshCellDecomp::debug)
    {
        Info<< "void tetPointMapperCellDecomp::calcAddressing() const : "
            << "Calculating addressing: ";
    }

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
        FatalErrorIn("void tetPointMapperCellDecomp::calcAddressing() const)")
            << "Addressing already calculated"
            << abort(FatalError);
    }

    const label oldCellOffset = mpm_.nOldPoints();

    // Mapping

    // Calculate direct (if all are direct)
    directPtr_ =
        new bool
        (
            pointMap_.direct()
         && cellMap_.direct()
        );

    // Assemble the maps
    if (*directPtr_)
    {
        // Direct mapping
        if (tetPolyMeshCellDecomp::debug)
        {
            Info<< " direct" << endl;
        }

        const labelList& mappedPoints = pointMap_.directAddressing();
        const labelList& mappedCells = cellMap_.directAddressing();

        directAddrPtr_ = new labelList(size());
        labelList& addr = *directAddrPtr_;
        label nAddr = 0;

        forAll (mappedPoints, pointI)
        {
            addr[nAddr] = mappedPoints[pointI];
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
        if (tetPolyMeshCellDecomp::debug)
        {
            Info<< " interpolative" << endl;
        }

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

            // Do cell addressing, 

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
         || cellMap_.insertedObjects()
        );

    // If there are, assemble the labels
    if (*insertedObjectsPtr_)
    {
        const labelList& insPoints = pointMap_.insertedObjectLabels();
        const labelList& insCells = cellMap_.insertedObjectLabels();

        insertedObjectLabelsPtr_ =
            new labelList(insPoints.size() + insCells.size());
        labelList& ins = *insertedObjectLabelsPtr_;

        label nIns = 0;

        forAll (insPoints, pointI)
        {
            ins[nIns] = insPoints[pointI];
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

    if (tetPolyMeshCellDecomp::debug)
    {
        Info<< "void tetPointMapperCellDecomp::calcAddressing() const : "
            << "Finished calculating addressing."
            << endl;
    }
}


void Foam::tetPointMapperCellDecomp::clearOut()
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
Foam::tetPointMapperCellDecomp::tetPointMapperCellDecomp
(
    const tetPolyMeshCellDecomp& mesh,
    const mapPolyMesh& meshMap,
    const pointMapper& pMapper,
    const cellMapper& cMapper
)
:
    mesh_(mesh),
    mpm_(meshMap),
    pointMap_(pMapper),
    cellMap_(cMapper),
    size_(mesh().nPoints() + mesh().nCells()),
    directPtr_(NULL),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    insertedObjectsPtr_(NULL),
    insertedObjectLabelsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tetPointMapperCellDecomp::~tetPointMapperCellDecomp()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::tetPointMapperCellDecomp::size() const
{
    return size_;
}


Foam::label Foam::tetPointMapperCellDecomp::sizeBeforeMapping() const
{
    return mpm_.nOldPoints() + mpm_.nOldCells();
}


bool Foam::tetPointMapperCellDecomp::direct() const
{
    if (!directPtr_)
    {
        calcAddressing();
    }

    return *directPtr_;
}


const Foam::unallocLabelList&
Foam::tetPointMapperCellDecomp::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const unallocLabelList& tetPointMapperCellDecomp::"
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


const Foam::labelListList& Foam::tetPointMapperCellDecomp::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& tetPointMapperCellDecomp::addressing() const"
        )   << "Requested interpolative addressing for a direct mapper."
            << abort(FatalError);
    }

    if (!interpolationAddrPtr_)
    {
        calcAddressing();
    }

    return *interpolationAddrPtr_;
}


const Foam::scalarListList& Foam::tetPointMapperCellDecomp::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& tetPointMapperCellDecomp::weights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcAddressing();
    }

    return *weightsPtr_;
}


bool Foam::tetPointMapperCellDecomp::insertedObjects() const
{
    if (!insertedObjectsPtr_)
    {
        calcAddressing();
    }

    return *insertedObjectsPtr_;
}


const Foam::labelList&
Foam::tetPointMapperCellDecomp::insertedObjectLabels() const
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
