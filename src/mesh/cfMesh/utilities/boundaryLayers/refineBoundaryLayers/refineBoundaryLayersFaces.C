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

#include "refineBoundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "demandDrivenData.H"
#include "FixedList.H"
#include "helperFunctions.H"

//#define DEBUGLayer

# ifdef DEBUGLayer
#include "OFstream.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void refineBoundaryLayers::refineFace
(
    const face& f,
    const FixedList<label, 2>& nLayersInDirection,
    DynList<DynList<label, 4>, 128>& newFaces
)
{
    //- this face must be a quad
    if( f.size() != 4 )
    {
        WarningIn
        (
            "void refineBoundaryLayers::refineFace(const face&,"
            " const FixedList<label, 2>&, DynList<DynList<label, 4> >&)"
        ) << "Face " << f << " is not a quad" << endl;
        return;
    }

    //- direction 0 represents edges 0 and 2
    //- direction 1 represents edges 1 and 3
    if( (nLayersInDirection[0] <= 1) && (nLayersInDirection[1] <= 1) )
    {
        //- this face may comprise of some split edges
        DynList<label, 64> newF;
        forAll(f, eI)
        {
            const edge e = f.faceEdge(eI);

            //- add the current point label
            newF.append(f[eI]);

            //- check if a split edge matches this face edge
            forAllRow(splitEdgesAtPoint_, f[eI], peI)
            {
                const label seI = splitEdgesAtPoint_(f[eI], peI);
                const edge& se = splitEdges_[seI];

                if( e == se )
                {
                    //- check the orientation and add new vertices created
                    //- on this edge
                    const label s = newVerticesForSplitEdge_.sizeOfRow(seI) - 1;
                    if( e.start() == se.start() )
                    {
                        for(label pI=1;pI<s;++pI)
                            newF.append(newVerticesForSplitEdge_(seI, pI));
                    }
                    else
                    {
                        for(label pI=s-1;pI>0;--pI)
                            newF.append(newVerticesForSplitEdge_(seI, pI));
                    }
                }
            }
        }

        newFaces.setSize(1);
        newFaces[0] = newF;
        return;
    }

    //- check which face edge is a direction 0 and which one is a direction 1
    label dir0(-1), dir1(-1);
    labelPair dir0Edges(-1, -1), dir1Edges(-1, -1);
    forAll(f, eI)
    {
        const edge e = f.faceEdge(eI);

        label ses(-1), see(-1);
        bool start(false), end(false);
        forAllRow(splitEdgesAtPoint_, e.start(), i)
        {
            const edge& se = splitEdges_[splitEdgesAtPoint_(e.start(), i)];

            if( (se.start() == e.start()) && (se.end() == f.prevLabel(eI)) )
            {
                ses = splitEdgesAtPoint_(e.start(), i);
                start = true;
                break;
            }
        }

        forAllRow(splitEdgesAtPoint_, e.end(), i)
        {
            const edge& se = splitEdges_[splitEdgesAtPoint_(e.end(), i)];

            if( (se.start() == e.end()) && (se.end() == f[(eI+2)%4]) )
            {
                see = splitEdgesAtPoint_(e.end(), i);
                end = true;
                break;
            }
        }

        if( start && end )
        {
            if( dir0 == -1 )
            {
                dir0 = eI;
                dir0Edges = labelPair(ses, see);
            }
            else if( dir1 == -1 )
            {
                dir1 = eI;
                dir1Edges = labelPair(ses, see);
            }
            else
            {
                FatalErrorIn
                (
                    "void refineBoundaryLayers::refineFace(const face&,"
                    " const FixedList<label, 2>&, DynList<DynList<label, 4> >&)"
                ) << "More than two split directions for a face"
                  << abort(FatalError);
            }
        }
    }

    # ifdef DEBUGLayer
    Pout << "Refining face " << f << endl;
    Pout << "Splits in direction " << nLayersInDirection << endl;
    Pout << "Here " << endl;
    Pout << "Dir0 " << dir0 << endl;
    Pout << "dir0Edges " << dir0Edges << endl;
    Pout << "Dir1 " << dir1 << endl;
    Pout << "dir1Edges " << dir1Edges << endl;
    # endif

    if( (dir0 < 0) && (dir1 < 0) )
    {
        Pout << "Refining face " << f << endl;
        forAll(f, pI)
        {
            if( splitEdgesAtPoint_.size() >= f[pI] )
            Pout << "Split edges at point " << f[pI]
                 << " are " << splitEdgesAtPoint_[f[pI]] << endl;
        }
        Pout << "Splits in direction " << nLayersInDirection << endl;
        Pout << "Here " << endl;
        Pout << "Dir0 " << dir0 << endl;
        Pout << "dir0Edges " << dir0Edges << endl;
        Pout << "Dir1 " << dir1 << endl;
        Pout << "dir1Edges " << dir1Edges << endl;

        FatalErrorIn
        (
            "void refineBoundaryLayers::refineFace(const face&,"
            " const FixedList<label, 2>&, DynList<DynList<label, 4> >&)"
        ) << "Cannot find split edges for a face" << abort(FatalError);
    }

    //- in case of only one refinement direction, it must direction 0
    if( (dir1 != -1) && (dir0 == -1) )
    {
        dir0 = dir1;
        dir0Edges = dir1Edges;
        dir1 = -1;
    }
    else if( (dir0 != -1) && (dir1 != -1) && (dir1 != f.fcIndex(dir0)) )
    {
        //- alternate value to preserve correct face orientation
        const label add = dir0;
        dir0 = dir1;
        dir1 = add;

        const labelPair lpAdd = dir0Edges;
        dir0Edges = dir1Edges;
        dir1Edges = lpAdd;
    }

    //- permutate the number of refinements in each direction
    const label nLayersDir0 = dir0>=0?nLayersInDirection[dir0%2]:1;
    const label nLayersDir1 = dir1>=0?nLayersInDirection[dir1%2]:1;

    # ifdef DEBUGLayer
    Pout << "Face has points " << f << endl;
    Pout << "dirEdges0 " << dir0Edges << endl;
    Pout << "dir1Edges " << dir1Edges << endl;
    if( dir0 >= 0 )
    {
        Pout << "Points on edge " << dir0Edges.first() << " with nodes "
             << splitEdges_[dir0Edges.first()]
             << " are " << newVerticesForSplitEdge_[dir0Edges.first()] << endl;
        Pout << "Points on edge " << dir0Edges.second() << " with nodes "
             << splitEdges_[dir0Edges.second()]
             << " are " << newVerticesForSplitEdge_[dir0Edges.second()] << endl;
    }
    if( dir1 >= 0 )
    {
        Pout << "Points on edge " << dir1Edges.first() << " with nodes "
             << splitEdges_[dir1Edges.first()]
             << " are " << newVerticesForSplitEdge_[dir1Edges.first()] << endl;
        Pout << "Points on edge " << dir1Edges.second() << " with nodes "
             << splitEdges_[dir1Edges.second()]
             << " are " << newVerticesForSplitEdge_[dir1Edges.second()] << endl;
    }
    Pout << "nLayersDir0 " << nLayersDir0 << endl;
    Pout << "nLayersDir1 " << nLayersDir1 << endl;
    # endif

    //- map the face onto a matrix for easier orientation
    DynList<DynList<label> > facePoints;
    facePoints.setSize(nLayersDir0+1);
    forAll(facePoints, i)
    {
        facePoints[i].setSize(nLayersDir1+1);
        facePoints[i] = -1;
    }

    //- add points in the matrix
    for(label i=0;i<nLayersDir0;++i)
    {
        facePoints[i][0] = newVerticesForSplitEdge_(dir0Edges.second(), i);
        facePoints[i][nLayersDir1] =
            newVerticesForSplitEdge_(dir0Edges.first(), i);
    }
    facePoints[nLayersDir0][0] = splitEdges_[dir0Edges.second()].end();
    facePoints[nLayersDir0][nLayersDir1] = splitEdges_[dir0Edges.first()].end();

    for(label i=1;i<nLayersDir1;++i)
    {
        facePoints[0][i] = newVerticesForSplitEdge_(dir1Edges.first(), i);
        facePoints[nLayersDir0][i] =
            newVerticesForSplitEdge_(dir1Edges.second(), i);
    }

    //- create missing vertices if there are any
    pointFieldPMG& points = mesh_.points();
    const point v00 = points[facePoints[0][0]];
    const point v10 = points[facePoints[nLayersDir0][0]];
    const point v01 = points[facePoints[0][nLayersDir1]];
    const point v11 = points[facePoints[nLayersDir0][nLayersDir1]];

    # ifdef DEBUGLayer
    forAll(points, pointI)
        Pout << "Point " << pointI << " coordinates " << points[pointI] << endl;
    Pout << "v00 = " << v00 << endl;
    Pout << "v10 = " << v10 << endl;
    Pout << "v11 = " << v11 << endl;
    Pout << "v01 = " << v01 << endl;
    # endif

    forAll(facePoints, i)
    {
        forAll(facePoints[i], j)
        {
            if( facePoints[i][j] < 0 )
            {
                # ifdef DEBUGLayer
                Pout << "Determining u " << facePoints[0][0]
                     << " coordinates " << points[facePoints[0][0]] << endl;
                Pout << "Other point " << facePoints[i][0]
                     << " coordinates " << points[facePoints[i][0]] << endl;
                Pout << "Points at aplit edge "
                     << newVerticesForSplitEdge_[dir0Edges.second()] << endl;
                # endif

                const scalar u
                (
                    Foam::mag(points[facePoints[i][0]] - v00) /
                    splitEdges_[dir0Edges.second()].mag(points)
                );

                # ifdef DEBUGLayer
                Pout << "Determining v " << facePoints[0][0]
                     << " coordinates " << points[facePoints[0][0]] << endl;
                Pout << "Other point " << facePoints[0][j]
                     << " coordinates " << points[facePoints[0][j]] << endl;
                Pout << "Points at aplit edge "
                     << newVerticesForSplitEdge_[dir1Edges.first()] << endl;
                # endif

                const scalar v
                (
                    Foam::mag(points[facePoints[0][j]] - v00) /
                    splitEdges_[dir1Edges.first()].mag(points)
                );

                # ifdef DEBUGLayer
                Pout << "Generating point of face " << endl;
                Pout << "u = " << u << endl;
                Pout << "v = " << v << endl;
                # endif

                //- calculate the coordinates of the missing point via
                //- transfinite interpolation
                const point newP
                (
                    (1.0 - u) * (1.0 - v) * v00 +
                    u * (1.0 - v) * v10 +
                    u * v * v11 +
                    (1.0 - u) * v * v01
                );

                # ifdef DEBUGLayer
                Pout << "Point coordinate " << newP << endl;
                # endif

                //- add the vertex to the mesh
                facePoints[i][j] = points.size();
                points.append(newP);
            }
        }
    }

    # ifdef DEBUGLayer
    Pout << "Face points after creating vertices " << facePoints << endl;
    # endif

    //- Finally, create the faces
    for(label j=0;j<nLayersDir1;++j)
    {
        for(label i=0;i<nLayersDir0;++i)
        {
            //- create quad face
            DynList<label, 4> f;

            f.append(facePoints[i][j]);

            if( (i == (nLayersDir0 - 1)) && (j == 0) )
            {
                # ifdef DEBUGLayer
                Pout << "1. Adding additional points on edge " << endl;
                # endif

                //- add additional points on edge
                const label eLabel = dir0Edges.second();
                const label size =
                    newVerticesForSplitEdge_.sizeOfRow(eLabel) - 1;

                for(label index=i+1;index<size;++index)
                    f.append(newVerticesForSplitEdge_(eLabel, index));
            }

            f.append(facePoints[i+1][j]);

            if(
                (dir1 != -1) &&
                (i == (nLayersDir0 - 1)) &&
                (j == (nLayersDir1 - 1))
            )
            {
                # ifdef DEBUGLayer
                Pout << "2. Adding additional points on edge " << endl;
                # endif

                //- add additional points on edge
                const label eLabel = dir1Edges.second();
                const label size =
                    newVerticesForSplitEdge_.sizeOfRow(eLabel) - 1;

                for(label index=j+1;index<size;++index)
                    f.append(newVerticesForSplitEdge_(eLabel, index));
            }

            f.append(facePoints[i+1][j+1]);

            if( (i == (nLayersDir0 - 1)) && (j == (nLayersDir1 - 1)) )
            {
                # ifdef DEBUGLayer
                Pout << "3. Adding additional points on edge " << endl;
                # endif

                const label eLabel = dir0Edges.first();
                const label size =
                    newVerticesForSplitEdge_.sizeOfRow(eLabel) - 2;
                for(label index=size;index>i;--index)
                    f.append(newVerticesForSplitEdge_(eLabel, index));
            }

            f.append(facePoints[i][j+1]);

            if( (dir1 != -1) && (i == 0) && (j == (nLayersDir1 - 1)) )
            {
                # ifdef DEBUGLayer
                Pout << "4. Adding additional points on edge " << endl;
                # endif

                const label eLabel = dir1Edges.first();
                const label size =
                    newVerticesForSplitEdge_.sizeOfRow(eLabel) - 2;
                for(label index=size;index>j;--index)
                    f.append(newVerticesForSplitEdge_(eLabel, index));
            }

            newFaces.append(f);
        }
    }

    # ifdef DEBUGLayer
    Pout << "Input face " << f << endl;
    Pout << "Decomposed faces are " << newFaces << endl;
    //if( (nLayersInDirection[0] > 1) && (nLayersInDirection[1] > 1) )
    //::exit(1);
    # endif
}

void refineBoundaryLayers::sortFacePoints
(
    const label faceI,
    DynList<DynList<label> >& facePoints,
    const label transpose
) const
{
    const faceListPMG& faces = mesh_.faces();
    const face& f = faces[faceI];

    # ifdef DEBUGLayer
    Pout << "Creating matrix of points on a split face " << faceI << endl;
    Pout << "Face comprises of points " << f << endl;
    Pout << "New faces from face " << facesFromFace_.sizeOfRow(faceI) << endl;
    # endif

    label procStart = mesh_.faces().size();
    const PtrList<processorBoundaryPatch>& procBoundaries = mesh_.procBoundaries();
    if( Pstream::parRun() )
        procStart = procBoundaries[0].patchStart();

    if(
        (faceI < procStart) ||
        procBoundaries[mesh_.faceIsInProcPatch(faceI)].owner()
    )
    {
        //- orientation of new faces is the same as the face itself
        //- start the procedure by finding the number of splits in
        //- both i and j direction
        label numSplitsI(1);

        const label pos = f.which(newFaces_(facesFromFace_(faceI, 0), 0));

        forAllRow(facesFromFace_, faceI, i)
        {
            const label nfI = facesFromFace_(faceI, i);

            if( (numSplitsI == 1) && newFaces_.contains(nfI, f.nextLabel(pos)) )
            {
                numSplitsI = i + 1;
                break;
            }
        }

        const label numSplitsJ = (facesFromFace_.sizeOfRow(faceI) / numSplitsI);

        # ifdef DEBUGLayer
        Pout << "Pos " << pos << endl;
        Pout << "Num splits in direction 0 " << numSplitsI << endl;
        Pout << "Num splits in direction 1 " << numSplitsJ << endl;
        # endif

        facePoints.setSize(numSplitsI+1);
        forAll(facePoints, i)
            facePoints[i].setSize(numSplitsJ+1);

        //- start filling in the matrix
        forAllRow(facesFromFace_, faceI, fI)
        {
            const label nfI = facesFromFace_(faceI, fI);

            const label i = fI % numSplitsI;
            const label j = fI / numSplitsI;

            # ifdef DEBUGLayer
            Pout << "New face " << fI << " is " << newFaces_[nfI] << endl;
            Pout << " i = " << i << endl;
            Pout << " j = " << j << endl;
            # endif

            if( newFaces_.sizeOfRow(nfI) == 4 )
            {
                facePoints[i][j] = newFaces_(nfI, 0);
                facePoints[i+1][j] = newFaces_(nfI, 1);
                facePoints[i+1][j+1] = newFaces_(nfI, 2);
                facePoints[i][j+1] = newFaces_(nfI, 3);
            }
            else
            {
                if( j == 0 )
                {
                    forAllRow(newFaces_, nfI, pI)
                        if( f.which(newFaces_(nfI, pI)) >= 0 )
                        {
                            facePoints[i+1][0] = newFaces_(nfI, pI);
                            break;
                        }
                }
                else if( i == 0 )
                {
                    forAllRow(newFaces_, nfI, pI)
                        if( f.which(newFaces_(nfI, pI)) >= 0 )
                        {
                            facePoints[0][j+1] = newFaces_(nfI, pI);
                            break;
                        }
                }
                else
                {
                    forAllRow(newFaces_, nfI, pI)
                        if( f.which(newFaces_(nfI, pI)) >= 0 )
                        {
                            facePoints[i+1][j+1] = newFaces_(nfI, pI);
                            break;
                        }
                }
            }
        }

        # ifdef DEBUGLayer
        Pout << "Generated matrix of points on face " << faceI
             << " is " << facePoints << endl;
        # endif
    }
    else
    {
        //- this situation exists on inter-processor boundaries
        //- on neighbour processor. i and j coordinates are reversed
        label numSplitsJ(1);

        const label pos = f.which(newFaces_(facesFromFace_(faceI, 0), 0));

        forAllRow(facesFromFace_, faceI, j)
        {
            const label nfI = facesFromFace_(faceI, j);

            if( (numSplitsJ == 1) && newFaces_.contains(nfI, f.prevLabel(pos)) )
            {
                numSplitsJ = j + 1;
                break;
            }
        }

        const label numSplitsI = (facesFromFace_.sizeOfRow(faceI) / numSplitsJ);

        # ifdef DEBUGLayer
        Pout << "2. Face comprises of points " << f << endl;
        Pout << "2. Num splits in direction 0 " << numSplitsI << endl;
        Pout << "2. Num splits in direction 1 " << numSplitsJ << endl;
        # endif

        facePoints.setSize(numSplitsI+1);
        forAll(facePoints, i)
            facePoints[i].setSize(numSplitsJ+1);

        //- start filling in the matrix
        forAllRow(facesFromFace_, faceI, fI)
        {
            const label nfI = facesFromFace_(faceI, fI);

            const label i = fI / numSplitsJ;
            const label j = fI % numSplitsJ;

            # ifdef DEBUGLayer
            Pout << "2. New face " << fI << " is " << newFaces_[nfI] << endl;
            Pout << "2. i = " << i << endl;
            Pout << "2. j = " << j << endl;
            # endif

            if( newFaces_.sizeOfRow(nfI) == 4 )
            {
                facePoints[i][j] = newFaces_(nfI, 0);
                facePoints[i+1][j] = newFaces_(nfI, 1);
                facePoints[i+1][j+1] = newFaces_(nfI, 2);
                facePoints[i][j+1] = newFaces_(nfI, 3);
            }
            else
            {
                if( i == 0 )
                {
                    forAllRow(newFaces_, nfI, pI)
                        if( f.which(newFaces_(nfI, pI)) >= 0 )
                        {
                            facePoints[0][j+1] = newFaces_(nfI, pI);
                            break;
                        }
                }
                else if( j == 0 )
                {
                    forAllRow(newFaces_, nfI, pI)
                        if( f.which(newFaces_(nfI, pI)) >= 0 )
                        {
                            facePoints[i+1][0] = newFaces_(nfI, pI);
                            break;
                        }
                }
                else
                {
                    forAllRow(newFaces_, nfI, pI)
                        if( f.which(newFaces_(nfI, pI)) >= 0 )
                        {
                            facePoints[i+1][j+1] = newFaces_(nfI, pI);
                            break;
                        }
                }
            }
        }

        # ifdef DEBUGLayer
        Pout << "Generated matrix of points on processor face " << faceI
             << " is " << facePoints << endl;
        # endif
    }

    if( transpose )
    {
        DynList<DynList<label> > transposedFacePoints;
        transposedFacePoints.setSize(facePoints[0].size());
        forAll(transposedFacePoints, j)
            transposedFacePoints[j].setSize(facePoints.size());

        forAll(facePoints, i)
            forAll(facePoints[i], j)
                transposedFacePoints[j][i] = facePoints[i][j];

        facePoints = transposedFacePoints;

        # ifdef DEBUGLayer
        Pout << "Transposed face points " << facePoints << endl;
        # endif
    }
}

void refineBoundaryLayers::sortFaceFaces
(
    const label faceI,
    DynList<DynList<label> >& faceFaces,
    const label transpose
) const
{
    const faceListPMG& faces = mesh_.faces();
    const face& f = faces[faceI];

    # ifdef DEBUGLayer
    Pout << "Creating matrix of faces on a split face " << faceI << endl;
    Pout << "Face comprises of points " << f << endl;
    # endif

    label procStart = mesh_.faces().size();
    const PtrList<processorBoundaryPatch>& procBoundaries = mesh_.procBoundaries();
    if( Pstream::parRun() )
        procStart = procBoundaries[0].patchStart();

    if(
        (faceI < procStart) ||
        procBoundaries[mesh_.faceIsInProcPatch(faceI)].owner()
    )
    {
        //- orientation of new faces is the same as the face itself
        //- start the procedure by finding the number of splits in
        //- both i and j direction
        label numSplitsI(1);

        const label pos = f.which(newFaces_(facesFromFace_(faceI, 0), 0));

        forAllRow(facesFromFace_, faceI, i)
        {
            const label nfI = facesFromFace_(faceI, i);

            if( (numSplitsI == 1) && newFaces_.contains(nfI, f.nextLabel(pos)) )
            {
                numSplitsI = i + 1;
                break;
            }
        }

        label numSplitsJ = (facesFromFace_.sizeOfRow(faceI) / numSplitsI);

        # ifdef DEBUGLayer
        Pout << "3. Num splits in direction 0 " << numSplitsI << endl;
        Pout << "3. Num splits in direction 1 " << numSplitsJ << endl;
        # endif

        faceFaces.setSize(numSplitsI);
        forAll(faceFaces, i)
            faceFaces[i].setSize(numSplitsJ);

        //- start filling in the matrix
        forAllRow(facesFromFace_, faceI, fI)
        {
            const label nfI = facesFromFace_(faceI, fI);

            const label i = fI % numSplitsI;
            const label j = fI / numSplitsI;

            # ifdef DEBUGLayer
            Pout << "3. New face " << fI << " is " << newFaces_[nfI] << endl;
            Pout << "3. i = " << i << endl;
            Pout << "3. j = " << j << endl;
            # endif

            faceFaces[i][j] = nfI;
        }

        # ifdef DEBUGLayer
        Pout << "3. Generated matrix of points on face " << faceI
             << " is " << faceFaces << endl;
        # endif
    }
    else
    {
        //- this situation exists on inter-processor boundaries
        //- on neighbour processor. i and j coordinates are reversed
        label numSplitsJ(1);

        const label pos = f.which(newFaces_(facesFromFace_(faceI, 0), 0));

        forAllRow(facesFromFace_, faceI, j)
        {
            const label nfI = facesFromFace_(faceI, j);

            if( (numSplitsJ == 1) && newFaces_.contains(nfI, f.prevLabel(pos)) )
            {
                numSplitsJ = j + 1;
                break;
            }
        }

        const label numSplitsI = (facesFromFace_.sizeOfRow(faceI) / numSplitsJ);

        # ifdef DEBUGLayer
        Pout << "4. Num splits in direction 0 " << numSplitsI << endl;
        Pout << "4. Num splits in direction 1 " << numSplitsJ << endl;
        # endif

        faceFaces.setSize(numSplitsI);
        forAll(faceFaces, i)
            faceFaces[i].setSize(numSplitsJ);

        //- start filling in the matrix
        forAllRow(facesFromFace_, faceI, fI)
        {
            const label nfI = facesFromFace_(faceI, fI);

            const label i = fI / numSplitsJ;
            const label j = fI % numSplitsJ;

            # ifdef DEBUGLayer
            Pout << "4. New face " << fI << " is " << newFaces_[nfI] << endl;
            Pout << "4. i = " << i << endl;
            Pout << "4. j = " << j << endl;
            # endif

            faceFaces[i][j] = nfI;
        }

        # ifdef DEBUGLayer
        Pout << "4. Generated matrix of faces on processor face " << faceI
             << " is " << faceFaces << endl;
        # endif
    }

    if( transpose )
    {
        DynList<DynList<label> > transposedFaceFaces;
        transposedFaceFaces.setSize(faceFaces[0].size());
        forAll(transposedFaceFaces, j)
            transposedFaceFaces[j].setSize(faceFaces.size());

        forAll(faceFaces, i)
            forAll(faceFaces[i], j)
                transposedFaceFaces[j][i] = faceFaces[i][j];

        faceFaces = transposedFaceFaces;

        # ifdef DEBUGLayer
        Pout << "Transposed face faces " << faceFaces << endl;
        # endif
    }
}

void refineBoundaryLayers::generateNewFaces()
{
    //- generate new boundary and inter-processor faces
    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& facePatches = mse.boundaryFacePatches();
    const edgeList& edges = mse.edges();
    const labelList& bp = mse.bp();
    const VRWGraph& bfEdges = mse.faceEdges();
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const VRWGraph& beFaces = mse.edgeFaces();

    //- mesh data
    const label nInternalFaces = mesh_.nInternalFaces();
    const faceListPMG& faces = mesh_.faces();

    //- container for faces
    facesFromFace_.setSize(faces.size());
    newFaces_.clear();

    //- split internal faces
    for(label faceI=0;faceI<nInternalFaces;++faceI)
    {
        const face& f = faces[faceI];

        //- only quad faces can be split
        if( f.size() != 4 )
        {
            facesFromFace_.append(faceI, newFaces_.size());
            newFaces_.appendList(f);
            continue;
        }

        //- check if there exist an edge of the face at the boundary
        FixedList<label, 2> nRefinementInDirection(1);

        forAll(f, eI)
        {
            const edge fe = f.faceEdge(eI);

            const label bps = bp[fe.start()];

            if( bps < 0 )
                continue;

            forAllRow(bpEdges, bps, bpsI)
            {
                const label beI = bpEdges(bps, bpsI);

                if( edges[beI] == fe )
                {
                    //- this edge is attached to the boundary
                    //- get the number of layers for neighbouring cells
                    const label nSplits0 = nLayersAtBndFace_[beFaces(beI, 0)];
                    const label nSplits1 = nLayersAtBndFace_[beFaces(beI, 1)];

                    //- set the number of layers for the given direction
                    const label dir = eI % 2;
                    nRefinementInDirection[dir] =
                        Foam::max
                        (
                            nRefinementInDirection[dir],
                            Foam::max(nSplits0, nSplits1)
                        );
                }
            }
        }

        //- refine the face
        DynList<DynList<label, 4>, 128> newFacesForFace;
        refineFace(f, nRefinementInDirection, newFacesForFace);

        //- store decomposed faces
        forAll(newFacesForFace, fI)
        {
            facesFromFace_.append(faceI, newFaces_.size());
            newFaces_.appendList(newFacesForFace[fI]);
        }

        # ifdef DEBUGLayer
        Pout << "Internal face " << faceI << " with points " << f
             << " is refined " << endl;
        forAllRow(facesFromFace_, faceI, i)
            Pout << "New face " << i << " is "
                 << newFaces_[facesFromFace_(faceI, i)] << endl;
        DynList<DynList<label> > tralala;
        sortFacePoints(faceI, tralala);
        # endif
    }

    //- refine boundary faces where needed
    //- it is required in locations where two or three layers intersect
    const label startingBoundaryFace = mesh_.boundaries()[0].patchStart();
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];
        const label faceI = startingBoundaryFace + bfI;

        //- only quad faces can be split
        if( bf.size() != 4 )
        {
            facesFromFace_.append(faceI, newFaces_.size());
            newFaces_.appendList(bf);
            continue;
        }

        //- check whether this face shall be refined and in which directions
        FixedList<label, 2> nRefinementInDirection(1);

        forAll(bf, eI)
        {
            const label beI = bfEdges(bfI, eI);

            if( beFaces.sizeOfRow(beI) != 2 )
                continue;

            //- get the neighbour face over the edge
            label neiFace = beFaces(beI, 0);

            if( neiFace == bfI )
                neiFace = beFaces(beI, 1);

            //- faces cannot be in the same layer
            const DynList<label>& neiLayers =
                layerAtPatch_[facePatches[neiFace]];

            if( neiLayers.size() == 0 )
                continue;

            const DynList<label>& currLayers = layerAtPatch_[facePatches[bfI]];

            bool foundSame(false);

            forAll(currLayers, i)
            {
                if( neiLayers.contains(currLayers[i]) )
                {
                    foundSame = true;
                    break;
                }
            }

            if( foundSame || (neiLayers.size() == 0) )
                continue;

            //- set the refinement direction for this face
            nRefinementInDirection[eI%2] = nLayersAtBndFace_[neiFace];
        }

        //- refine the face
        DynList<DynList<label, 4>, 128> newFacesForFace;
        refineFace(bf, nRefinementInDirection, newFacesForFace);

        //- store the refined faces
        forAll(newFacesForFace, fI)
        {
            facesFromFace_.append(faceI, newFaces_.size());
            newFaces_.appendList(newFacesForFace[fI]);
        }

        # ifdef DEBUGLayer
        Pout << "Boundary face " << faceI << " with points " << bf
             << " owner cell " << mesh_.owner()[faceI] << " is refined " << endl;
        forAllRow(facesFromFace_, faceI, i)
            Pout << "New face " << i << " is "
                 << newFaces_[facesFromFace_(faceI, i)] << endl;
        # endif
    }

    if( Pstream::parRun() )
    {
        //- refine faces at interprocessor boundaries
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();

        //- exchange information about the number of splits
        //- to other processors
        std::map<label, DynList<labelPair, 2> > localSplits;
        forAll(procBoundaries, patchI)
        {
            labelLongList sendData;

            const label start = procBoundaries[patchI].patchStart();
            const label size = procBoundaries[patchI].patchSize();

            for(label fI=0;fI<size;++fI)
            {
                const label faceI = start + fI;
                const face& f = faces[faceI];

                forAll(f, eI)
                {
                    const edge fe = f.faceEdge(eI);

                    const label bps = bp[fe.start()];

                    if( bps < 0 )
                        continue;

                    forAllRow(bpEdges, bps, bpeI)
                    {
                        const label beI = bpEdges(bps, bpeI);

                        if( edges[beI] == fe )
                        {
                            //- this edge is attached to the boundary
                            //- get the number of layers for neighbouring cell
                            const label nSplits0 =
                                nLayersAtBndFace_[beFaces(beI, 0)];

                            //- add the data to the list for sending
                            const label dir = (eI % 2);

                            # ifdef DEBUGLayer
                            Pout << "Face " << fI << " owner of proc patch "
                                 << procBoundaries[patchI].myProcNo()
                                 << " nei proc "
                                 << procBoundaries[patchI].neiProcNo()
                                 << " bnd face patch "
                                 << facePatches[beFaces(beI, 0)]
                                 << " direction " << dir
                                 << " nSplits " << nSplits0 << endl;
                            # endif

                            //- add face label, direction
                            //- and the number of splits
                            sendData.append(fI);
                            sendData.append(dir);
                            sendData.append(nSplits0);
                            localSplits[faceI].append(labelPair(dir, nSplits0));
                        }
                    }
                }
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                sendData.byteSize()
            );

            toOtherProc << sendData;
        }

        //- receive data from other procesors
        forAll(procBoundaries, patchI)
        {
            //- get the data sent from the neighbour processor
            labelList receivedData;

            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            fromOtherProc >> receivedData;

            const label start = procBoundaries[patchI].patchStart();

            label counter(0);
            while( counter < receivedData.size() )
            {
                const label fI = receivedData[counter++];
                const label dir = ((receivedData[counter++] + 1) % 2);
                const label nSplits = receivedData[counter++];

                DynList<labelPair, 2>& currentSplits = localSplits[start+fI];
                forAll(currentSplits, i)
                {
                    if( currentSplits[i].first() == dir )
                        currentSplits[i].second() =
                            Foam::max(currentSplits[i].second(), nSplits);
                }
            }
        }

        # ifdef DEBUGLayer
        returnReduce(1, sumOp<label>());
        Pout << "Starting splitting processor boundaries" << endl;
        # endif

        //- perform splitting
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label size = procBoundaries[patchI].patchSize();

            for(label fI=0;fI<size;++fI)
            {
                const label faceI = start + fI;

                std::map<label, DynList<labelPair, 2> >::const_iterator it =
                    localSplits.find(faceI);

                if( it == localSplits.end() )
                {
                    //- this face is not split
                    facesFromFace_.append(faceI, newFaces_.size());
                    newFaces_.appendList(faces[faceI]);
                    continue;
                }

                //- split the face and add the faces to the list
                if( procBoundaries[patchI].owner() )
                {
                    //- this processor owns this patch
                    FixedList<label, 2> nLayersInDirection(1);
                    const DynList<labelPair, 2>& dirSplits = it->second;
                    forAll(dirSplits, i)
                        nLayersInDirection[dirSplits[i].first()] =
                            dirSplits[i].second();

                    # ifdef DEBUGLayer
                    Pout << "Face " << fI << " at owner processor "
                        << procBoundaries[patchI].myProcNo()
                        << " neighbour processor "
                        << procBoundaries[patchI].neiProcNo()
                        << " face " << faces[faceI] << " refinement direction "
                        << nLayersInDirection << endl;
                    # endif

                    DynList<DynList<label, 4>, 128> facesFromFace;
                    refineFace(faces[faceI], nLayersInDirection, facesFromFace);

                    //- add faces
                    forAll(facesFromFace, i)
                    {
                        facesFromFace_.append(faceI, newFaces_.size());
                        newFaces_.appendList(facesFromFace[i]);
                    }
                }
                else
                {
                    //- reverse the face before splitting
                    FixedList<label, 2> nLayersInDirection(1);
                    const DynList<labelPair, 2>& dirSplits = it->second;
                    forAll(dirSplits, i)
                        nLayersInDirection[(dirSplits[i].first()+1)%2] =
                            dirSplits[i].second();

                    const face rFace = faces[faceI].reverseFace();

                    # ifdef DEBUGLayer
                    Pout << "Face " << fI << " at owner processor "
                        << procBoundaries[patchI].myProcNo()
                        << " neighbour processor "
                        << procBoundaries[patchI].neiProcNo()
                        << " face " << rFace << " refinement direction "
                        << nLayersInDirection << endl;
                    # endif

                    DynList<DynList<label, 4>, 128> facesFromFace;
                    refineFace(rFace, nLayersInDirection, facesFromFace);

                    forAll(facesFromFace, i)
                    {
                        const DynList<label, 4>& df = facesFromFace[i];
                        DynList<label, 4> rFace = help::reverseFace(df);

                        facesFromFace_.append(faceI, newFaces_.size());
                        newFaces_.appendList(rFace);
                    }
                }
            }
        }

        # ifdef DEBUGLayer
        returnReduce(1, sumOp<label>());
        for(label procI=0;procI<Pstream::nProcs();++procI)
        {
            if( procI == Pstream::myProcNo() )
            {
                forAll(procBoundaries, patchI)
                {
                    const label start = procBoundaries[patchI].patchStart();
                    const label size = procBoundaries[patchI].patchSize();

                    for(label fI=0;fI<size;++fI)
                    {
                        const label faceI = start + fI;
                        const face& f = faces[faceI];
                        Pout << "Face " << fI << " in patch "
                             << procBoundaries[patchI].patchName()
                             << " has nodes " << f
                             << " local splits " << localSplits[faceI]
                             << " new faces from face " << facesFromFace_[faceI]
                             << endl;

                        Pout << " Face points ";
                        forAll(f, pI)
                            Pout << mesh_.points()[f[pI]] << " ";
                        Pout << endl;

                        forAllRow(facesFromFace_, faceI, ffI)
                        {
                            const label nfI = facesFromFace_(faceI, ffI);
                            Pout << "New face " << ffI << " with label " << nfI
                                 << " consists of points ";
                            forAllRow(newFaces_, nfI, pI)
                                Pout << mesh_.points()[newFaces_(nfI, pI)]
                                     << " ";
                            Pout << endl;
                        }
                    }
                }
            }

            returnReduce(1, sumOp<label>());
        }

        returnReduce(1, sumOp<label>());
        //::exit(1);
        # endif
    }

    # ifdef DEBUGLayer
    returnReduce(1, sumOp<label>());

    for(label procI=0;procI<Pstream::nProcs();++procI)
    {
        if( procI == Pstream::myProcNo() )
        {
            Pout << "facesFromFace_ " << facesFromFace_ << endl;
            Pout << "newFaces_ " << newFaces_ << endl;
        }

        returnReduce(1, sumOp<label>());
    }

    OFstream file("refinedFaces.vtk");

    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    //- write points
    file << "POINTS " << mesh_.points().size() << " float\n";
    forAll(mesh_.points(), pI)
    {
        const point& p = mesh_.points()[pI];

        file << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
    }

    //- write faces
    label counter(0);
    forAll(newFaces_, faceI)
    {
        counter += newFaces_.sizeOfRow(faceI);
        ++counter;
    }

    file << "\nPOLYGONS " << faces.size()
         << " " << counter << nl;
    forAll(newFaces_, faceI)
    {
        file << newFaces_.sizeOfRow(faceI);
        forAllRow(newFaces_, faceI, i)
            file << " " << newFaces_(faceI, i);
        file << nl;
    }

    file << "\n";
    # endif

    Info << "Finished refining boundary-layer faces " << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
