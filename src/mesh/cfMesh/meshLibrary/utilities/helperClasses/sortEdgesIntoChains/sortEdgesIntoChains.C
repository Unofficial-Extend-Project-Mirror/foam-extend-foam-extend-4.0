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

#include "sortEdgesIntoChains.H"
#include "helperFunctions.H"
#include "Map.H"

//#define DEBUGSort

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void sortEdgesIntoChains::createNodeLabels()
{
    label nPoints(0);
    forAll(bEdges_, eI)
    {
        const edge& e = bEdges_[eI];
        if( !newNodeLabel_.found(e.start()) )
            newNodeLabel_.insert(e.start(), nPoints++);
        if( !newNodeLabel_.found(e.end()) )
            newNodeLabel_.insert(e.end(), nPoints++);
    }

    edgesAtPoint_.setSize(nPoints, DynList<label>());
    forAll(bEdges_, eI)
    {
        const edge& e = bEdges_[eI];
        label l = newNodeLabel_[e.start()];
        edgesAtPoint_[l].append(eI);

        l = newNodeLabel_[e.end()];
        edgesAtPoint_[l].append(eI);
    }

    forAll(edgesAtPoint_, pI)
        if( edgesAtPoint_[pI].size() % 2 )
            openEdges_ = true;
}

bool sortEdgesIntoChains::findPointsBelongingToTheChain
(
    const label currPos,
    boolList& chainEdges
) const
{
    # ifdef DEBUGSort
    Info << "Finding point belonging to a chain" << endl;
    # endif

    chainEdges.setSize(bEdges_.size());
    chainEdges = false;

    if( edgesAtPoint_[currPos].size() != 2 )
        return false;

    const label commonVrt =
        bEdges_[edgesAtPoint_[currPos][0]].commonVertex
            (
                bEdges_[edgesAtPoint_[currPos][1]]
            );
    label prevVrt = bEdges_[edgesAtPoint_[currPos][0]].otherVertex(commonVrt);
    label nextVrt = bEdges_[edgesAtPoint_[currPos][1]].otherVertex(commonVrt);
    forAll(edgesAtPoint_[currPos], posI)
        chainEdges[edgesAtPoint_[currPos][posI]] = true;

    # ifdef DEBUGSort
    Info << "commonVrt " << commonVrt << endl;
    Info << "prevVrt " << prevVrt << endl;
    Info << "nextVrt " << nextVrt << endl;
    # endif

    bool found;
    do
    {
        found = false;

        const DynList<label>& vEdges = edgesAtPoint_[newNodeLabel_[prevVrt]];
        if( vEdges.size() == 2 )
        {
            forAll(vEdges, eI)
                if( !chainEdges[vEdges[eI]] )
                {
                    found = true;
                    chainEdges[vEdges[eI]] = true;
                    prevVrt = bEdges_[vEdges[eI]].otherVertex(prevVrt);
                }
        }
    } while( found );

    do
    {
        found = false;

        const DynList<label>& vEdges = edgesAtPoint_[newNodeLabel_[nextVrt]];
        if( vEdges.size() == 2 )
        {
            forAll(vEdges, eI)
                if( !chainEdges[vEdges[eI]] )
                {
                    found = true;
                    chainEdges[vEdges[eI]] = true;
                    nextVrt = bEdges_[vEdges[eI]].otherVertex(nextVrt);
                }
        }
    } while( found );

    if(
        (edgesAtPoint_[newNodeLabel_[nextVrt]].size() != 2) &&
        (edgesAtPoint_[newNodeLabel_[prevVrt]].size() != 2) &&
        (prevVrt != nextVrt)
    )
    {
        chainEdges = false;
        return false;
    }

    # ifdef DEBUGSort
    Info << "Chain edges " << chainEdges << endl;
    # endif

    return true;
}

void sortEdgesIntoChains::shrinkEdges(const boolList& chainEdges)
{
    forAll(chainEdges, eI)
        if( chainEdges[eI] )
        {
            const edge& e = bEdges_[eI];
            edgesAtPoint_[newNodeLabel_[e.start()]].removeElement
            (
                edgesAtPoint_[newNodeLabel_[e.start()]].containsAtPosition(eI)
            );

            edgesAtPoint_[newNodeLabel_[e.end()]].removeElement
            (
                edgesAtPoint_[newNodeLabel_[e.end()]].containsAtPosition(eI)
            );
        }
}

void sortEdgesIntoChains::createChainFromEdges(const boolList& chainEdges)
{
    direction i(0);
    forAll(chainEdges, eI)
        if( chainEdges[eI] )
            ++i;

    labelList chainPoints(i);
    i = 0;

    forAll(chainEdges, eI)
        if( chainEdges[eI] )
        {
            chainPoints[i++] = bEdges_[eI].start();
            chainPoints[i++] = bEdges_[eI].end();

            # ifdef DEBUGSort
            Info << "Init chainPoints " << chainPoints << endl;
            # endif

            bool found;
            do
            {
                # ifdef DEBUGSort
                Info << "Iteration " << label(i-1) << endl;
                # endif

                found = false;
                const DynList<label>& pEdges =
                    edgesAtPoint_[newNodeLabel_[chainPoints[i-1]]];

                forAll(pEdges, peI)
                    if( chainEdges[pEdges[peI]] )
                    {
                        const label otherPoint =
                            bEdges_[pEdges[peI]].otherVertex(chainPoints[i-1]);

                        # ifdef DEBUGSort
                        Info << "Other point " << otherPoint << endl;
                        # endif
                        if( otherPoint == -1 )
                            continue;
                        if( chainPoints[i-2] == otherPoint )
                            continue;
                        if( chainPoints[0] == otherPoint )
                            continue;

                        found = true;
                        chainPoints[i++] = otherPoint;
                    }
            } while( found );

            createdChains_.append(chainPoints);

            break;
        }
}

void sortEdgesIntoChains::sortEdges()
{
    createNodeLabels();

    if( !openEdges_ )
    {
        boolList chainEdges(bEdges_.size());
        forAll(edgesAtPoint_, pI)
            if( findPointsBelongingToTheChain(pI, chainEdges) )
            {
                createChainFromEdges(chainEdges);

                shrinkEdges(chainEdges);
            }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

sortEdgesIntoChains::sortEdgesIntoChains(const DynList<edge>& bEdges)
:
    bEdges_(bEdges),
    openEdges_(false),
    newNodeLabel_(),
    edgesAtPoint_(),
    createdChains_()
{
    sortEdges();
}

sortEdgesIntoChains::~sortEdgesIntoChains()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
// Member functions
const DynList<labelList>& sortEdgesIntoChains::sortedChains() const
{
    return createdChains_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //
