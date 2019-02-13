/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "meshOctreeModifier.H"
#include "HashSet.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define OCTREETiming
//#define OCTREE_DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeModifier::updateCommunicationPattern()
{
    # ifdef OCTREETiming
    const scalar startTime = omp_get_wtime();
    # endif

    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    //- create the list which contains ranges of addresses at a given processor
    List<Pair<meshOctreeCubeCoordinates> > range(Pstream::nProcs());

    //- create the range for the current processor
    range[Pstream::myProcNo()].first() = leaves[0]->coordinates();
    range[Pstream::myProcNo()].second() =
        leaves[leaves.size()-1]->coordinates();

    //- communicate missing ranges
    Pstream::gatherList(range);
    Pstream::scatterList(range);

    //- find missing child cubes in the tree. These coordinates are located on
    //- other processors, and they must fit in the range of cubes located
    //- on that processor.
    LongList<meshOctreeCubeCoordinates> migratedCubes;
    octree_.initialCubePtr_->findCoordinatesOfMissingCubes(migratedCubes);

    labelHashSet newNeiProcs;
    forAll(migratedCubes, mcI)
    {
        const meshOctreeCubeCoordinates& cc = migratedCubes[mcI];
        forAll(range, procI)
        {
            if( procI == Pstream::myProcNo() )
                continue;

            if( (cc >= range[procI].first()) && (cc <= range[procI].second()) )
                newNeiProcs.insert(procI);
        }
    }

    //- create new neighbour procs and their range
    octree_.neiProcs_ = newNeiProcs.toc();

    octree_.neiRange_.setSize(octree_.neiProcs_.size());
    forAll(octree_.neiProcs_, i)
        octree_.neiRange_[i] = range[octree_.neiProcs_[i]];

    # ifdef OCTREETiming
    Info << "Time for updating communication pattern "
        << (omp_get_wtime()-startTime) << endl;
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
