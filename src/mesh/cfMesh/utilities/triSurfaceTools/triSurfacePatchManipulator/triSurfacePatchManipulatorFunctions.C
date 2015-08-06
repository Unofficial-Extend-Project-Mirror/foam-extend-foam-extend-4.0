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

#include "triSurfacePatchManipulator.H"
#include "helperFunctions.H"
#include "demandDrivenData.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfacePatchManipulator::allocateFeatureEdges()
{
    const edgeLongList& edges = surf_.edges();
    const VRWGraph& pEdges = surf_.pointEdges();

    //- allocate featureEdges list
    featureEdges_.setSize(edges.size());
    featureEdges_ = direction(0);

    const edgeLongList& featureEdges = surf_.featureEdges();

    forAll(featureEdges, feI)
    {
        const edge& e = featureEdges[feI];

        forAllRow(pEdges, e.start(), peI)
        {
            const label eI = pEdges(e.start(), peI);

            if( edges[eI] == e )
                featureEdges_[eI] |= 1;
        }
    }
}

void triSurfacePatchManipulator::createPatches()
{
    nPatches_ = 0;
    facetInPatch_.setSize(surf_.size());
    facetInPatch_ = -1;

    const VRWGraph& faceEdges = surf_.facetEdges();
    const VRWGraph& edgeFaces = surf_.edgeFacets();

    forAll(facetInPatch_, triI)
    {
        if( facetInPatch_[triI] != -1 )
            continue;

        labelLongList front;
        front.append(triI);
        facetInPatch_[triI] = nPatches_;

        while( front.size() )
        {
            const label fLabel = front.removeLastElement();

            const constRow fEdges = faceEdges[fLabel];

            forAll(fEdges, feI)
            {
                const label edgeI = fEdges[feI];

                //- check if th edges is marked as a feature edge
                if( featureEdges_[edgeI] )
                    continue;

                const constRow eFaces = edgeFaces[edgeI];

                //- stop at non-manifold edges
                if( eFaces.size() != 2 )
                    continue;

                label neiTri = eFaces[0];
                if( neiTri == fLabel )
                    neiTri = eFaces[1];

                //- do not overwrite existing patch information
                if( surf_[fLabel].region() != surf_[neiTri].region() )
                    continue;
                if( facetInPatch_[neiTri] != -1 )
                    continue;

                facetInPatch_[neiTri] = nPatches_;
                front.append(neiTri);
            }
        }

        ++nPatches_;
    }

    Info << "Created " << nPatches_ << " surface patches" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
