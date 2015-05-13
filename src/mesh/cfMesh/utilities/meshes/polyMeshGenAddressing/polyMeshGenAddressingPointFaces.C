/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"
#include "VRWGraphSMPModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcPointFaces() const
{
    if( pfPtr_)
    {
        FatalErrorIn("polyMeshGenAddressing::calcPointFaces()")
            << "pointFaces already calculated"
            << abort(FatalError);
    }
    else
    {
        const faceListPMG& faces = mesh_.faces();
        const label nPoints = mesh_.points().size();

        //- create the addressing
        pfPtr_ = new VRWGraph();
        VRWGraph& pointFaceAddr = *pfPtr_;

        VRWGraphSMPModifier(pointFaceAddr).reverseAddressing(faces);
        pointFaceAddr.setSize(nPoints);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const VRWGraph& polyMeshGenAddressing::pointFaces() const
{
    if( !pfPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& polyMeshGenAddressing::pointFaces() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcPointFaces();
    }

    return *pfPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
