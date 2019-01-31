/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "checkCellConnectionsOverFaces.H"
#include "labelLongList.H"
#include "labelledPair.H"

# ifdef USE_OMP
#include <omp.h>
# endif

#include <set>
#include <map>

#include "helperFunctions.H"

//#define DEBUGCheck

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private member functions * * * * * * * * * * * * * * * //

namespace meshConnectionsHelper
{

class meshConnectionsNeighbourOperator
{
    const polyMeshGen& mesh_;

public:

    meshConnectionsNeighbourOperator(const polyMeshGen& mesh)
    :
        mesh_(mesh)
    {}

    label size() const
    {
        return mesh_.cells().size();
    }

    void operator()(const label cellI, DynList<label>& neighbourCells) const
    {
        neighbourCells.clear();

        const labelList& owner = mesh_.owner();
        const labelList& neighbour = mesh_.neighbour();

        const cell& c = mesh_.cells()[cellI];

        forAll(c, fI)
        {
            label nei = owner[c[fI]];

            if( nei == cellI )
                nei = neighbour[c[fI]];

            if( nei >= 0 )
                neighbourCells.append(nei);
        }
    }

    template<class labelListType>
    void collectGroups
    (
        std::map<label, DynList<label> >& neiGroups,
        const labelListType& elementInGroup,
        const DynList<label>& localGroupLabel
    ) const
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();
        const labelList& owner = mesh_.owner();

        //- send the data to other processors
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label size = procBoundaries[patchI].patchSize();

            labelList groupOwner(procBoundaries[patchI].patchSize());
            for(label faceI=0;faceI<size;++faceI)
            {
                const label groupI = elementInGroup[owner[start+faceI]];

                if( groupI < 0 )
                {
                    groupOwner[faceI] = -1;
                    continue;
                }

                groupOwner[faceI] = localGroupLabel[groupI];
            }

            OPstream toOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo(),
                groupOwner.byteSize()
            );

            toOtherProc << groupOwner;
        }

        //- receive data from other processors
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();

            labelList receivedData;

            IPstream fromOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            fromOtherProc >> receivedData;

            forAll(receivedData, faceI)
            {
                if( receivedData[faceI] < 0 )
                    continue;

                const label groupI = elementInGroup[owner[start+faceI]];

                if( groupI < 0 )
                    continue;

                DynList<label>& ng = neiGroups[localGroupLabel[groupI]];

                //- store the connection over the inter-processor boundary
                ng.appendIfNotIn(receivedData[faceI]);
            }
        }
    }
};

class meshConnectionsSelectorOperator
{

public:

    meshConnectionsSelectorOperator()
    {}

    bool operator()(const label /*cellI*/) const
    {
        return true;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meshConnectionsHelper

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkCellConnectionsOverFaces::findCellGroups()
{
    Info << "Checking cell connections" << endl;

    mesh_.owner();
    nGroups_ =
        help::groupMarking
        (
            cellGroup_,
            meshConnectionsHelper::meshConnectionsNeighbourOperator(mesh_),
            meshConnectionsHelper::meshConnectionsSelectorOperator()
        );

    Info << "Finished checking cell connections" << endl;
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //
// Constructors

checkCellConnectionsOverFaces::checkCellConnectionsOverFaces(polyMeshGen& mesh)
:
    mesh_(mesh),
    cellGroup_(mesh.cells().size(), -1),
    nGroups_(0)
{
    findCellGroups();
}

// * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * * * * //

checkCellConnectionsOverFaces::~checkCellConnectionsOverFaces()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool checkCellConnectionsOverFaces::checkCellGroups()
{
    if( nGroups_ == 1 )
        return false;

    Warning << "Mesh has " << nGroups_ << " unconnected regions" << endl;

    labelList nCellsInGroup(nGroups_, 0);

    forAll(cellGroup_, cI)
        ++nCellsInGroup[cellGroup_[cI]];

    if( Pstream::parRun() )
    {
        forAll(nCellsInGroup, groupI)
            reduce(nCellsInGroup[groupI], sumOp<label>());
    }

    //- find groups which has most cells this group will be kept
    label maxGroup(-1);
    forAll(nCellsInGroup, groupI)
        if( nCellsInGroup[groupI] > maxGroup )
        {
            maxGroup = nCellsInGroup[groupI];
            nGroups_ = groupI;
        }

    //- remove cells which are not in the group which has max num of cells
    boolList removeCell(mesh_.cells().size(), false);
    forAll(cellGroup_, cellI)
        if( cellGroup_[cellI] != nGroups_ )
            removeCell[cellI] = true;

    polyMeshGenModifier(mesh_).removeCells(removeCell);

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
