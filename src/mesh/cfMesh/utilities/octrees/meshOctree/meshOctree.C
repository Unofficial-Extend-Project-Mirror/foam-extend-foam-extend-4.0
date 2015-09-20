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
#include "triSurf.H"
#include "boundBox.H"
#include "demandDrivenData.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from surface
meshOctree::meshOctree(const triSurf& ts, const bool isQuadtree)
:
    surface_(ts),
    neiProcs_(),
    neiRange_(),
    initialCubePtr_(NULL),
    initialCubeRotation_(0),
    rootBox_(),
    isRootInitialised_(false),
    searchRange_(0.0),
    octantVectors_(),
    vrtLeavesPos_(),
    regularityPositions_(),
    dataSlots_(),
    leaves_(),
    isQuadtree_(isQuadtree)
{
    createInitialOctreeBox();

    setOctantVectorsAndPositions();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshOctree::~meshOctree()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctree::setOctantVectorsAndPositions()
{
    octantVectors_[0] = Vector<label>(-1, -1, -1);
    octantVectors_[1] = Vector<label>(1, -1, -1);
    octantVectors_[2] = Vector<label>(-1, 1, -1);
    octantVectors_[3] = Vector<label>(1, 1, -1);
    octantVectors_[4] = Vector<label>(-1, -1, 1);
    octantVectors_[5] = Vector<label>(1, -1, 1);
    octantVectors_[6] = Vector<label>(-1, 1, 1);
    octantVectors_[7] = Vector<label>(1, 1, 1);

    # ifdef DEBUGSearch
    Info << "Octant vectors " << octantVectors_ << endl;
    # endif

    //- set regularity positions
    //- neighbours over faces
    regularityPositions_[0] = meshOctreeCubeCoordinates(-1, 0, 0, 0);
    regularityPositions_[1] = meshOctreeCubeCoordinates(1, 0, 0, 0);
    regularityPositions_[2] = meshOctreeCubeCoordinates(0, -1, 0, 0);
    regularityPositions_[3] = meshOctreeCubeCoordinates(0, 1, 0, 0);
    regularityPositions_[4] = meshOctreeCubeCoordinates(0, 0, -1, 0);
    regularityPositions_[5] = meshOctreeCubeCoordinates(0, 0, 1, 0);

    //- neighbours over edges
    //- edges in x-direction
    regularityPositions_[6] = meshOctreeCubeCoordinates(0, -1, -1, 0);
    regularityPositions_[7] = meshOctreeCubeCoordinates(0, 1, -1, 0);
    regularityPositions_[8] = meshOctreeCubeCoordinates(0, -1, 1, 0);
    regularityPositions_[9] = meshOctreeCubeCoordinates(0, 1, 1, 0);

    //- edges in y-direction
    regularityPositions_[10] = meshOctreeCubeCoordinates(-1, 0, -1, 0);
    regularityPositions_[11] = meshOctreeCubeCoordinates(1, 0, -1, 0);
    regularityPositions_[12] = meshOctreeCubeCoordinates(-1, 0, 1, 0);
    regularityPositions_[13] = meshOctreeCubeCoordinates(1, 0, 1, 0);

    //- edges in z-direction
    regularityPositions_[14] = meshOctreeCubeCoordinates(-1, -1, 0, 0);
    regularityPositions_[15] = meshOctreeCubeCoordinates(1, -1, 0, 0);
    regularityPositions_[16] = meshOctreeCubeCoordinates(-1, 1, 0, 0);
    regularityPositions_[17] = meshOctreeCubeCoordinates(1, 1, 0, 0);

    //- neighbours over vertices
    regularityPositions_[18] = meshOctreeCubeCoordinates(-1, -1, -1, 0);
    regularityPositions_[19] = meshOctreeCubeCoordinates(1, -1, -1, 0);
    regularityPositions_[20] = meshOctreeCubeCoordinates(-1, 1, -1, 0);
    regularityPositions_[21] = meshOctreeCubeCoordinates(1, 1, -1, 0);
    regularityPositions_[22] = meshOctreeCubeCoordinates(-1, -1, 1, 0);
    regularityPositions_[23] = meshOctreeCubeCoordinates(1, -1, 1, 0);
    regularityPositions_[24] = meshOctreeCubeCoordinates(-1, 1, 1, 0);
    regularityPositions_[25] = meshOctreeCubeCoordinates(1, 1, 1, 0);

    # ifdef DEBUGSearch
    Info << "Regularity positions " << regularityPositions_ << endl;
    # endif

    //- set vrtLeavesPos_
    for(label vrtI=0;vrtI<8;++vrtI)
    {
        FixedList<label, 3> vc(0);

        if( vrtI & 1 )
            vc[0] += 1;
        if( vrtI & 2 )
            vc[1] += 1;
        if( vrtI & 4 )
            vc[2] += 1;

        # ifdef DEBUGSearch
        Info << "Vert " << vrtI << " vc " << vc << endl;
        # endif

        for(label i=0;i<8;++i)
        {
            FixedList<label, 3> pos;

            for(label j=0;j<3;++j)
            {
                if( vc[j] == 0 && octantVectors_[i][j] == 1 )
                {
                    pos[j] = 0;
                }
                else if( vc[j] == 1 && octantVectors_[i][j] == 1 )
                {
                    pos[j] = 1;
                }
                else
                {
                    pos[j] = vc[j] + octantVectors_[i][j];
                }
            }

            vrtLeavesPos_[vrtI][i] =
                meshOctreeCubeCoordinates
                (
                    pos[0],
                    pos[1],
                    pos[2],
                    direction(0)
                );
        }
    }

    # ifdef DEBUGSearch
    Info << "vrtLeavesPos_ " << vrtLeavesPos_ << endl;
    # endif
}

void meshOctree::createInitialOctreeBox()
{
    //- create initial octree box
    boundBox bb(surface_.points());
    const point& min_ = bb.min();
    const point& max_ = bb.max();

    const point c = (max_ + min_) / 2.0;
    scalar cs = 1.5 * (max_.x() - min_.x()) / 2.0;
    if( cs < (1.5 * (max_.y() - min_.y()) / 2.0) )
    {
        cs = 1.5 * (max_.y() - min_.y()) / 2.0;
    }
    if( cs < (1.5 * (max_.z() - min_.z()) / 2.0) )
    {
        cs = 1.5 * (max_.z() - min_.z()) / 2.0;
    }

    //- create root box and initial cube
    rootBox_ = boundBox(c - point(cs, cs, cs), c + point(cs, cs, cs));

    if( Pstream::parRun() )
    {
        reduce(rootBox_.min(), minOp<point>());
        reduce(rootBox_.max(), maxOp<point>());
    }

    //- allocate data slots
    # ifdef USE_OMP
    if( omp_get_num_procs() > 0 )
    {
        dataSlots_.setSize(omp_get_num_procs());
    }
    else
    {
        dataSlots_.setSize(1);
    }
    # else
    dataSlots_.setSize(1);
    # endif

    meshOctreeSlot* slotPtr = &dataSlots_[0];

    if( !isQuadtree_ )
    {
        slotPtr->cubes_.append
        (
            meshOctreeCube
            (
                meshOctreeCubeCoordinates(0, 0, 0, 0),
                surface_.size(),
                slotPtr
            )
        );
    }
    else
    {
        slotPtr->cubes_.append
        (
            meshOctreeCube
            (
                meshOctreeCubeCoordinates(0, 0, -10, 0),
                surface_.size(),
                slotPtr
            )
        );
    }

    initialCubePtr_ = &slotPtr->cubes_[0];

    leaves_.setSize(1);
    leaves_[0] = initialCubePtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
