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

#include "meshOctreeModifier.H"

//#define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeModifier::distributeLeavesToProcessors()
{
    if( octree_.neiProcs().size() != 0 )
        return;

    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;
    const label nLeavesPerProcessor = leaves.size() / Pstream::nProcs();

    Info << "Number of leaves per processor " << nLeavesPerProcessor << endl;

    if( nLeavesPerProcessor == 0 )
        return;

    Info << "Distributing leaves to processors" << endl;
    //- leaf boxes are sorted in Z-ordering
    //- initial distribution is performed by decomposing the list of leaves
    //- into the equal-size chunks and allocate them to processors
    List<Pair<meshOctreeCubeCoordinates> > procRange(Pstream::nProcs());
    label balance = leaves.size() % Pstream::nProcs();
    label start = leaves.size() - 1;
    for(short procI=Pstream::nProcs()-1;procI>=0;--procI)
    {
        label end = start - nLeavesPerProcessor;
        if( balance )
        {
            --balance;
            --end;
        }
        else if( procI == 0 )
        {
            end = -1;
        }

        procRange[procI].second() = leaves[start]->coordinates();
        procRange[procI].first() = leaves[end+1]->coordinates();

        for(;start>end;--start)
            leaves[start]->setProcNo(procI);
    }

    # ifdef DEBUGSearch
    forAll(leaves, leafI)
        Info << "Leaf " << leafI << " is in proc "
            << leaves[leafI]->procNo() << endl;
    # endif

    //- find neighbouring processors of the current processor
    DynList<label> neiProcs;

    forAll(leaves, leafI)
    {
        const meshOctreeCubeCoordinates& cc = leaves[leafI]->coordinates();

        forAll(procRange, procI)
        {
            if( procI == Pstream::myProcNo() )
                continue;

            if(
                (cc >= procRange[procI].first()) &&
                (cc <= procRange[procI].second())
            )
            {
                neiProcs.appendIfNotIn(procI);
                break;
            }
        }
    }

    # ifdef DEBUGSearch
    Serr << "Neighbouring processors of the processor " << Pstream::myProcNo()
        << " are " << neiProcs << endl;
    # endif

    labelList& nProcs = octree_.neiProcs_;
    nProcs.setSize(neiProcs.size());
    forAll(nProcs, i)
        nProcs[i] = neiProcs[i];
    Foam::sort(nProcs);

    octree_.neiRange_.transfer(procRange);

    # ifdef DEBUGSearch
    Info << "Nei procs " << octree_.neiProcs_ << endl;
    Info << "Nei range " << octree_.neiRange_ << endl;
    # endif

    //- delete cubes which are not local to the processor
    octree_.initialCubePtr_->purgeProcessorCubes(Pstream::myProcNo());

    //- update the list of leaves
    createListOfLeaves();

    Info << "Finished distributing leaves to processors" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
