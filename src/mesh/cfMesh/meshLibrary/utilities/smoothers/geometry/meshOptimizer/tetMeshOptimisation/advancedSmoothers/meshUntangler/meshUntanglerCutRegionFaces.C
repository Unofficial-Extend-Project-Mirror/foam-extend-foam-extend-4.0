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

#include "demandDrivenData.H"
#include "meshUntangler.H"

//#define DEBUGSmooth

# ifdef DEBUGSmooth
#include "helperFunctions.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshUntangler::cutRegion::findNewFaces()
{
    #ifdef DEBUGSmooth
    Info << "Finding new faces " << endl;
    #endif

    bool finished;
    do
    {
        finished = true;

        const DynList<DynList<label, 8>, 64>& fcs = *facesPtr_;
        DynList<edge, 128>& cEdges = *cEdgesPtr_;

        cFacesPtr_ = new DynList<DynList<label, 8>, 64>();
        DynList<DynList<label, 8>, 64>& cFaces = *cFacesPtr_;

        DynList<label, 8> faceInPlane;

        DynList<label, 64> pointUsage;
        pointUsage.setSize(cPtsPtr_->size());

        forAll(fcs, fI)
        {
            const DynList<label, 8>& f = fcs[fI];

            # ifdef DEBUGSmooth
            Info << "Creating new face from face " << fI
                << " consisting of edges " << f << endl;
            # endif

            pointUsage = 0;

            DynList<label, 8> newFace;

            forAll(f, eI)
            {
                # ifdef DEBUGSmooth
                const DynList<edge>& edges = *edgesPtr_;
                Info << "Vertex types for face edge " << eI << " are "
                    << label(vertexTypes_[edges[f[eI]].start()]) << " and "
                    << label(vertexTypes_[edges[f[eI]].end()]) << endl;
                # endif

                const label edgeLabel = newEdgeLabel_[f[eI]];

                if( edgeLabel != -1 )
                {
                    # ifdef DEBUGSmooth
                    Info << "Orig edge " << eI << " " << edges[f[eI]]
                        << " is replaced with " << cEdges[edgeLabel] << endl;
                    # endif

                    const edge& e = cEdges[edgeLabel];
                    ++pointUsage[e[0]];
                    ++pointUsage[e[1]];
                    newFace.append(edgeLabel);
                }
            }

            if( newFace.size() > 1 )
            {
                DynList<label, 4> newEdge;
                forAll(pointUsage, pI)
                    if( pointUsage[pI] == 1 )
                        newEdge.append(pI);

                if( newEdge.size() == 2 )
                {
                    # ifdef DEBUGSmooth
                    Info << "Storing new edge " << newEdge << endl;
                    # endif

                    newFace.append(cEdges.size());
                    cEdges.append(edge(newEdge[0], newEdge[1]));
                }
                else if( newEdge.size() > 2 )
                {
                    # ifdef DEBUGSmooth
                    Info << "New edge " << newEdge << endl;
                    # endif

                    tieBreak(f);
                    if( !valid_ ) return;
                    finished = false;
                    break;

                    FatalErrorIn
                    (
                        "void meshUntangler::cutRegion::findNewFaces()"
                    ) << "Edge has more than two nodes!"
                        << abort(FatalError);
                }

                cFaces.append(newFace);
            }
        }

        if( !finished ) continue;

        //- find edges which form the faceInPlane
        DynList<label, 128> edgeUsage;
        edgeUsage.setSize(cEdges.size());
        edgeUsage = 0;
        forAll(cFaces, fI)
        {
            const DynList<label, 8>& f = cFaces[fI];

            forAll(f, eI)
                ++edgeUsage[f[eI]];
        }

        forAll(edgeUsage, eI)
            if( edgeUsage[eI] == 1 )
                faceInPlane.append(eI);

        if( faceInPlane.size() > 2 )
        {
            # ifdef DEBUGSmooth
            Info << "Adding face in plane " << faceInPlane << endl;
            Info << "Face in plane consists of edges " << endl;
            forAll(faceInPlane, eI)
                Info << "Edge " << eI << " is "
                    << cEdges[faceInPlane[eI]] << endl;
            # endif

            cFaces.append(faceInPlane);
        }

        # ifdef DEBUGSmooth
        Info << "cEdges " << cEdges << endl;
        Info << "Number of faces before cutting " << fcs.size() << endl;
        Info << "Found " << cFaces.size() << " new faces" << endl;
        forAll(fcs, fI)
        {
            Info << "Old face " << fI << " contains edges " << fcs[fI] << endl;
        }

        forAll(cFaces, fI)
        {
            Info << "New face " << fI << " contains edges " << cFaces[fI] << endl;
        }

        //- test if the region is closed
        List<DynList<label, 4> > eFaces(cEdges.size());
        forAll(cFaces, fI)
        {
            const DynList<label, 8>& f = cFaces[fI];
            forAll(f, eI)
                eFaces[f[eI]].append(fI);
        }

        if( eFaces.size() > 5 )
        forAll(eFaces, fI)
            if( eFaces[fI].size() != 2 )
            {
                Info << "eFaces " << eFaces << endl;
                Info << "cEdges " << cEdges << endl;
                Info << "cFaces " << cFaces << endl;

                FatalErrorIn
                (
                    "void meshOptimizer::meshUntangler::"
                    "cutRegion::findNewFaces()"
                ) << "Cell is not topologically closed!" << abort(FatalError);
            }
        # endif

    } while( !finished );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
