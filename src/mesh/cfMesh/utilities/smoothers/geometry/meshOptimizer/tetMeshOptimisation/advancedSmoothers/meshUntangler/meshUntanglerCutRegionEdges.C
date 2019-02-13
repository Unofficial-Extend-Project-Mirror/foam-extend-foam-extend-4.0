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

#include "demandDrivenData.H"
#include "meshUntangler.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshUntangler::cutRegion::findNewEdges()
{
    #ifdef DEBUGSmooth
    Info << "Creating new edgesPtr_" << endl;
    #endif

    cEdgesPtr_ = new DynList<edge, 128>();
    DynList<edge, 128>& cEdges = *cEdgesPtr_;

    const DynList<point, 64>& pts = *pointsPtr_;
    const DynList<edge, 128>& edges = *edgesPtr_;

    newEdgeLabel_.setSize(edges.size());
    newEdgeLabel_ = label(-1);

    forAll(edges, eI)
    {
        const edge& e = edges[eI];
        const label start = e.start();
        const label end = e.end();
        const label newStartLabel = newVertexLabel_[start];
        const label newEndLabel = newVertexLabel_[end];

        if( (newStartLabel != -1) && (newEndLabel != -1) )
        {
            newEdgeLabel_[eI] = cEdges.size();
            cEdges.append(edge(newStartLabel, newEndLabel));
        }
        else if( (newEndLabel != -1) && (vertexTypes_[end] & KEEP) )
        {
            //- start edge vertex is not visible, but the other one is
            newEdgeLabel_[eI] = cEdges.size();
            cEdges.append(edge(newEndLabel, cPtsPtr_->size()));

            const scalar t =
                -vertexDistance_[start] /
                (vertexDistance_[end] - vertexDistance_[start]);

            const point newP = (1.0 - t) * pts[start] + t * pts[end];
            cPtsPtr_->append(newP);
        }
        else if( (newStartLabel != -1) && (vertexTypes_[start] & KEEP) )
        {
            //- end edge vertex is not visible, but the other one is
            newEdgeLabel_[eI] = cEdges.size();
            cEdges.append(edge(newStartLabel, cPtsPtr_->size()));

            const scalar t =
                -vertexDistance_[end] /
                (vertexDistance_[start] - vertexDistance_[end]);

            const point newP = (1.0 - t) * pts[end] + t * pts[start];
            cPtsPtr_->append(newP);
        }

        # ifdef DEBUGSmooth
        if( newEdgeLabel_[eI] != -1 )
        {
            Info << nl << "Edge " << eI << "consisting of " << e
                << " has been replaced with edge " << newEdgeLabel_[eI] << "  "
                << cEdges[newEdgeLabel_[eI]] << endl;
        }
        else
        {
            Info << "Edge " << e << " has been deleted!" << endl;
        }
        # endif
    }

    #ifdef DEBUGSmooth
    Info << "Found " << cEdges.size() << " new edges" << endl;
    #endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
