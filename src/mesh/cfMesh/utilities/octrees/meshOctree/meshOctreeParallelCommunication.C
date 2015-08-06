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

#include "meshOctree.H"
#include "HashSet.H"

//#define OCTREE_DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctree::exchangeRequestsWithNeighbourProcessors
(
    const LongList<meshOctreeCubeCoordinates>& dataToSend,
    LongList<meshOctreeCubeCoordinates>& dataToReceive
) const
{
    if( !Pstream::parRun() || (neiProcs_.size() == 0) )
    {
        //- this is a serial run
        return;
    }

    List<LongList<meshOctreeCubeCoordinates> > toProcs
    (
        neiProcs_.size(),
        LongList<meshOctreeCubeCoordinates>()
    );
    meshOctreeCubeCoordinates minCoord, maxCoord;
    forAll(dataToSend, i)
    {
        dataToSend[i].neighbourRange(minCoord, maxCoord);

        # ifdef OCTREE_DEBUG
        label counter(0);
        # endif

        forAll(neiProcs_, procI)
        {
            if( maxCoord >= neiRange_[procI].first() )
            {
                if( minCoord <= neiRange_[procI].second() )
                {
                    toProcs[procI].append(dataToSend[i]);

                    # ifdef OCTREE_DEBUG
                    ++counter;
                    # endif
                }
            }
            else
            {
                break;
            }
        }

        # ifdef OCTREE_DEBUG
        if( !counter )
        {
            DynList<label> neighs;
            this->findAllLeafNeighbours(dataToSend[i], neighs);
            forAll(neighs, nI)
            if( neighs[nI] == meshOctreeCubeBasic::OTHERPROC )
            {
                forAll(neiRange_, j)
                {
                    if( minCoord <= neiRange_[j].second() )
                        Pout << "Min coord " << minCoord << " is smaller "
                            << "than " << neiRange_[j].first() << endl;
                    if( maxCoord >= neiRange_[j].first() )
                        Pout << "Max coord " << maxCoord << " is greater "
                            << "than " << neiRange_[j].second() << endl;
                }

                FatalError << "Box " << dataToSend[i] << " is not sent"
                    << " minCoord " << minCoord
                    << " maxCoord " << maxCoord
                    << " neighbours " << neiProcs_
                    << " neiRange_ " << neiRange_
                    << abort(FatalError);
            }
        }
        # endif
    }

    //- send an receive the size of data chunk which will be exchanged
    forAll(neiProcs_, neiProcI)
    {
        OPstream toOtherProc
        (
            Pstream::blocking,
            neiProcs_[neiProcI],
            sizeof(label)
        );

        toOtherProc << toProcs[neiProcI].size();
    }

    labelList sizeOfOtherProc(neiProcs_.size());
    forAll(neiProcs_, neiProcI)
    {
        IPstream fromOtherProc
        (
            Pstream::blocking,
            neiProcs_[neiProcI],
            sizeof(label)
        );

        fromOtherProc >> sizeOfOtherProc[neiProcI];
    }

    //- exchange data between processors
    //- upper-diagonal communication
    forAll(neiProcs_, neiProcI)
    {
        if( sizeOfOtherProc[neiProcI] == 0 )
            continue;
        if( neiProcs_[neiProcI] >= Pstream::myProcNo() )
            continue;

        //- receive data from other processor
        IPstream fromOtherProc(Pstream::scheduled, neiProcs_[neiProcI]);
        dataToReceive.appendFromStream(fromOtherProc);
    }

    forAll(neiProcs_, neiProcI)
    {
        if( toProcs[neiProcI].size() == 0 )
            continue;
        if( neiProcs_[neiProcI] <= Pstream::myProcNo() )
            continue;

        //- send data to other processor
        OPstream toOtherProc(Pstream::scheduled, neiProcs_[neiProcI]);
        toOtherProc << toProcs[neiProcI];
    }

    //- lower-diagonal communication
    forAllReverse(neiProcs_, neiProcI)
    {
        if( sizeOfOtherProc[neiProcI] == 0 )
            continue;
        if( neiProcs_[neiProcI] <= Pstream::myProcNo() )
            continue;

        //- receive data from other processor
        IPstream fromOtherProc(Pstream::scheduled, neiProcs_[neiProcI]);
        dataToReceive.appendFromStream(fromOtherProc);
    }

    forAllReverse(neiProcs_, neiProcI)
    {
        if( toProcs[neiProcI].size() == 0 )
            continue;
        if( neiProcs_[neiProcI] >= Pstream::myProcNo() )
            continue;

        //- send data to other processor
        OPstream toOtherProc(Pstream::scheduled, neiProcs_[neiProcI]);
        toOtherProc << toProcs[neiProcI];
    }
}
void meshOctree::exchangeRequestsWithNeighbourProcessors
(
    const LongList<meshOctreeCubeCoordinates>& dataToSend,
    const LongList<scalar>& rangesToSend,
    LongList<meshOctreeCubeCoordinates>& dataToReceive,
    LongList<scalar>& receivedRanges
) const
{
    if( !Pstream::parRun() || (neiProcs_.size() == 0) )
    {
        //- this is a serial run
        return;
    }

    List<LongList<meshOctreeCubeCoordinates> > toProcs
    (
        neiProcs_.size(),
        LongList<meshOctreeCubeCoordinates>()
    );
    List<LongList<scalar> > attributesToProcs
    (
        neiProcs_.size(),
        LongList<scalar>()
    );

    meshOctreeCubeCoordinates minCoord, maxCoord;
    forAll(dataToSend, i)
    {
        dataToSend[i].neighbourRange(minCoord, maxCoord);

        const scalar size = dataToSend[i].size(rootBox_);
        const label nLayers = ceil(rangesToSend[i] / size);

        minCoord =
            meshOctreeCubeCoordinates
            (
                minCoord.posX() - nLayers,
                minCoord.posY() - nLayers,
                minCoord.posZ() - nLayers,
                minCoord.level()
            );
        maxCoord =
            meshOctreeCubeCoordinates
            (
                nLayers + maxCoord.posX(),
                nLayers + maxCoord.posY(),
                nLayers + maxCoord.posZ(),
                maxCoord.level()
            );

        # ifdef OCTREE_DEBUG
        label counter(0);
        # endif

        forAll(neiProcs_, procI)
        {
            if( maxCoord >= neiRange_[procI].first() )
            {
                if( minCoord <= neiRange_[procI].second() )
                {
                    toProcs[procI].append(dataToSend[i]);
                    attributesToProcs[procI].append(rangesToSend[i]);

                    # ifdef OCTREE_DEBUG
                    ++counter;
                    # endif
                }
            }
            else
            {
                break;
            }
        }
    }

    //- send an receive the size of data chunk which will be exchanged
    forAll(neiProcs_, neiProcI)
    {
        OPstream toOtherProc
        (
            Pstream::blocking,
            neiProcs_[neiProcI],
            sizeof(label)
        );

        toOtherProc << toProcs[neiProcI].size();
    }

    labelList sizeOfOtherProc(neiProcs_.size());
    forAll(neiProcs_, neiProcI)
    {
        IPstream fromOtherProc
        (
            Pstream::blocking,
            neiProcs_[neiProcI],
            sizeof(label)
        );

        fromOtherProc >> sizeOfOtherProc[neiProcI];
    }

    //- exchange data between processors
    //- upper-diagonal communication
    forAll(neiProcs_, neiProcI)
    {
        if( sizeOfOtherProc[neiProcI] == 0 )
            continue;
        if( neiProcs_[neiProcI] >= Pstream::myProcNo() )
            continue;

        //- receive data from other processor
        IPstream fromOtherProc(Pstream::scheduled, neiProcs_[neiProcI]);

        dataToReceive.appendFromStream(fromOtherProc);
    }

    forAll(neiProcs_, neiProcI)
    {
        if( toProcs[neiProcI].size() == 0 )
            continue;
        if( neiProcs_[neiProcI] <= Pstream::myProcNo() )
            continue;

        //- send data to other processor
        OPstream toOtherProc(Pstream::scheduled, neiProcs_[neiProcI]);
        toOtherProc << toProcs[neiProcI];
    }

    //- lower-diagonal communication
    forAllReverse(neiProcs_, neiProcI)
    {
        if( sizeOfOtherProc[neiProcI] == 0 )
            continue;
        if( neiProcs_[neiProcI] <= Pstream::myProcNo() )
            continue;

        //- receive data from other processor
        IPstream fromOtherProc(Pstream::scheduled, neiProcs_[neiProcI]);

        dataToReceive.appendFromStream(fromOtherProc);
    }

    forAllReverse(neiProcs_, neiProcI)
    {
        if( toProcs[neiProcI].size() == 0 )
            continue;
        if( neiProcs_[neiProcI] >= Pstream::myProcNo() )
            continue;

        //- send data to other processor
        OPstream toOtherProc(Pstream::scheduled, neiProcs_[neiProcI]);
        toOtherProc << toProcs[neiProcI];
    }

    //- exchange attributes
    //- exchange data between processors
    //- upper-diagonal communication
    forAll(neiProcs_, neiProcI)
    {
        if( sizeOfOtherProc[neiProcI] == 0 )
            continue;
        if( neiProcs_[neiProcI] >= Pstream::myProcNo() )
            continue;

        //- receive data from other processor
        IPstream fromOtherProc(Pstream::scheduled, neiProcs_[neiProcI]);

        receivedRanges.appendFromStream(fromOtherProc);
    }

    forAll(neiProcs_, neiProcI)
    {
        if( toProcs[neiProcI].size() == 0 )
            continue;
        if( neiProcs_[neiProcI] <= Pstream::myProcNo() )
            continue;

        //- send data to other processor
        OPstream toOtherProc
        (
            Pstream::scheduled,
            neiProcs_[neiProcI],
            attributesToProcs[neiProcI].byteSize()
        );

        toOtherProc << attributesToProcs[neiProcI];
    }

    //- lower-diagonal communication
    forAllReverse(neiProcs_, neiProcI)
    {
        if( sizeOfOtherProc[neiProcI] == 0 )
            continue;
        if( neiProcs_[neiProcI] <= Pstream::myProcNo() )
            continue;

        //- receive data from other processor
        IPstream fromOtherProc(Pstream::scheduled, neiProcs_[neiProcI]);

        receivedRanges.appendFromStream(fromOtherProc);
    }

    forAllReverse(neiProcs_, neiProcI)
    {
        if( toProcs[neiProcI].size() == 0 )
            continue;
        if( neiProcs_[neiProcI] >= Pstream::myProcNo() )
            continue;

        //- send data to other processor
        OPstream toOtherProc
        (
            Pstream::scheduled,
            neiProcs_[neiProcI],
            attributesToProcs[neiProcI].byteSize()
        );

        toOtherProc << attributesToProcs[neiProcI];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
