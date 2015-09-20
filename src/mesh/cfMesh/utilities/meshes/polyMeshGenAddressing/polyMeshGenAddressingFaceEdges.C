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

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcFaceEdges() const
{
    if( fePtr_ )
    {
        FatalErrorIn("polyMeshGenAddressing::calcFaceEdges() const")
            << "faceEdges already calculated"
            << abort(FatalError);
    }
    else
    {
        fePtr_ = new VRWGraph();
        VRWGraph& faceEdgesAddr = *fePtr_;

        const edgeList& edges = this->edges();

        const VRWGraph& pointFaces = this->pointFaces();
        const faceListPMG& faces = mesh_.faces();

        labelList nfe(faces.size());

        # ifdef USE_OMP
        const label nThreads = 3 * omp_get_num_procs();
        # endif

        # ifdef USE_OMP
        # pragma omp parallel num_threads(nThreads) if( faces.size() > 10000 )
        # endif
        {
            # ifdef USE_OMP
            # pragma omp for schedule(static)
            # endif
            forAll(faces, faceI)
                nfe[faceI] = faces[faceI].size();

            # ifdef USE_OMP
            # pragma omp barrier

            # pragma omp master
            # endif
            VRWGraphSMPModifier(faceEdgesAddr).setSizeAndRowSize(nfe);

            # ifdef USE_OMP
            # pragma omp barrier

            # pragma omp for schedule(static)
            # endif
            forAll(edges, edgeI)
            {
                const edge ee = edges[edgeI];
                const label s = ee.start();

                forAllRow(pointFaces, s, pfI)
                {
                    const label faceI = pointFaces(s, pfI);

                    const face& f = faces[faceI];
                    forAll(f, eI)
                    {
                        if( f.faceEdge(eI) == ee )
                        {
                            faceEdgesAddr[faceI][eI] = edgeI;
                            break;
                        }
                    }
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const VRWGraph& polyMeshGenAddressing::faceEdges() const
{
    if( !fePtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& polyMeshGenAddressing::faceEdges() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcFaceEdges();
    }

    return *fePtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
