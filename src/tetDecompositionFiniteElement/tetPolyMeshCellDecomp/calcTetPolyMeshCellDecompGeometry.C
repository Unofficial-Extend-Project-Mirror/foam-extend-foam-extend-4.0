/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "tetPolyMeshCellDecomp.H"
#include "tetPointRef.H"
#include "SquareMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label tetPolyMeshCellDecomp::nPoints() const
{
    if (nPoints_ < 0)
    {
        nPoints_ = mesh_.nPoints() + mesh_.nCells();
    }

    return nPoints_;
}


label tetPolyMeshCellDecomp::nEdges() const
{
    if (nEdges_ < 0)
    {
        nEdges_ = mesh_.nEdges();

        // add edges going across the face
        const faceList& f  = mesh_.faces();

        forAll (f, fI)
        {
            nEdges_ += f[fI].size() - 3;
        }

        // add point-to-cell edges
        const labelListList& pc = mesh_.pointCells();

        forAll (pc, pointI)
        {
            nEdges_ += pc[pointI].size();
        }
    }

    return nEdges_;
}


label tetPolyMeshCellDecomp::nTets() const
{
    if (nTets_ < 0)
    {
        // Count the cells
        nTets_ = 0;

        const cellList& polyCells = mesh_.cells();

        forAll (polyCells, cellI)
        {
            nTets_ += nTetsForCell(cellI);
        }
    }

    return nTets_;
}


tmp<pointField> tetPolyMeshCellDecomp::points() const
{
    tmp<pointField> ttetPoints(new pointField(nPoints()));
    pointField& tetPoints = ttetPoints();

    const pointField& points = mesh_.points();

    const pointField& cellCentres = mesh_.cellCentres();

    label tetPointI = 0;

    forAll (points, pointI)
    {
        tetPoints[tetPointI] = points[pointI];
        tetPointI++;
    }

    forAll (cellCentres, cellI)
    {
        tetPoints[tetPointI] = cellCentres[cellI];
        tetPointI++;
    }

    return ttetPoints;
}


cellShapeList tetPolyMeshCellDecomp::tetCells() const
{
    cellShapeList t(nTets());

    const cellList& polyCells = mesh_.cells();

    label nTetCells = 0;

    forAll (polyCells, cellI)
    {
        tetCellList cellTets = tets(cellI);

        forAll (cellTets, tetI)
        {
            t[nTetCells] = cellTets[tetI].tetCellShape();
            nTetCells++;
        }
    }

    return t;
}


// Return tetrahedral decomposition for cell
tetCellList tetPolyMeshCellDecomp::tets(const label cellID) const
{
    const unallocFaceList& meshFaces = mesh_.faces();

    const unallocLabelList& owner = mesh_.faceOwner();

    // Initialise the size of the return
    tetCellList t(nTetsForCell(cellID));

    label nTetras = 0;

    const labelList& cellFaces = mesh_.cells()[cellID];

    forAll (cellFaces, faceI)
    {
        const label curFaceID = cellFaces[faceI];

        face f = meshFaces[curFaceID];

        // Take care of owner/neighbour face reversal
        if (cellID == owner[curFaceID])
        {
            f = f.reverseFace();
        }

        // Grab the point zero of the face and create triangles by grabbing the
        // labels around the perimeter
        label pointZeroID = f[0];

        for (label i = 1; i < f.size() - 1; i++)
        {
            // Face is assumed to be inward-pointing
            t[nTetras] =
                tetCell
                (
                    pointZeroID,
                    f[i],
                    f[i + 1],
                    cellID + cellOffset_
                );

            nTetras++;
        }
    }

    return t;
}


void tetPolyMeshCellDecomp::gradNiDotGradNj
(
    const label cellID,
    SquareMatrix<scalar>& denseMatrix,
    const labelList& globalToLocalBuffer
) const
{
    const unallocFaceList& meshFaces = mesh_.faces();
    const unallocLabelList& owner = mesh_.faceOwner();
    const labelList& cellFaces = mesh_.cells()[cellID];

    const pointField& points = mesh_.points();
    const pointField& cellCentres = mesh_.cellCentres();

    scalarField tetBuffer(tetPointRef::nEdges);

    // For optimisation, point allocation is taken out of the loop
    label vertexI, edgeI, localI, localJ;
    edge curEdge;

    forAll (cellFaces, faceI)
    {
        const label curFaceID = cellFaces[faceI];

        face f = meshFaces[curFaceID];

        if (cellID == owner[curFaceID])
        {
            f = f.reverseFace();
        }

        label pointZeroID = f[0];

        for (label i = 1; i < f.size() - 1; i++)
        {
            tetCell curTet
            (
                pointZeroID,
                f[i],
                f[i + 1],
                cellOffset() + cellID
            );

            tetPointRef tpr
            (
                points[pointZeroID],
                points[f[i]],
                points[f[i + 1]],
                cellCentres[cellID]
            );

            // Calculate the diagonal and insert into local matrix
            tpr.gradNiSquared(tetBuffer);

            for
            (
                vertexI = 0;
                vertexI < tetPointRef::nVertices;
                vertexI++
            )
            {
                localI = globalToLocalBuffer[curTet[vertexI]];

                denseMatrix[localI][localI] += tetBuffer[vertexI];
            }

            // Note: the coefficient is symmetric, i.e. equal in the upper and
            // lower triangle
            tpr.gradNiDotGradNj(tetBuffer);

            for (edgeI = 0; edgeI < tetPointRef::nEdges; edgeI++)
            {
                edge curEdge = curTet.tetEdge(edgeI);

                localI = globalToLocalBuffer[curEdge.start()];
                localJ = globalToLocalBuffer[curEdge.end()];

                // Insert the coefficient into upper and lower triangle
                // HJ, 4/Dec/2006
                denseMatrix[localI][localJ] += tetBuffer[edgeI];
                denseMatrix[localJ][localI] += tetBuffer[edgeI];
            }
        }
    }
}


void tetPolyMeshCellDecomp::gradNiGradNj
(
    const label cellID,
    SquareMatrix<tensor>& denseMatrix,
    const labelList& globalToLocalBuffer
) const
{
    const unallocFaceList& meshFaces = mesh_.faces();
    const unallocLabelList& owner = mesh_.faceOwner();
    const labelList& cellFaces = mesh_.cells()[cellID];

    const pointField& points = mesh_.points();
    const pointField& cellCentres = mesh_.cellCentres();

    tensorField tetBuffer(tetPointRef::nEdges);

    // For optimisation, point allocation is taken out of the loop
    label vertexI, edgeI, localI, localJ, minIJ, maxIJ;
    edge curEdge;

    forAll (cellFaces, faceI)
    {
        const label curFaceID = cellFaces[faceI];

        face f = meshFaces[curFaceID];

        if (cellID == owner[curFaceID])
        {
            f = f.reverseFace();
        }

        label pointZeroID = f[0];

        for (label i = 1; i < f.size() - 1; i++)
        {
            tetCell curTet
            (
                pointZeroID,
                f[i],
                f[i + 1],
                cellOffset() + cellID
            );

            tetPointRef tpr
            (
                points[pointZeroID],
                points[f[i]],
                points[f[i + 1]],
                cellCentres[cellID]
            );

            // Calculate the diagonal and insert into local matrix
            tpr.gradNiGradNi(tetBuffer);

            for
            (
                vertexI = 0;
                vertexI < tetPointRef::nVertices;
                vertexI++
            )
            {
                localI = globalToLocalBuffer[curTet[vertexI]];

                denseMatrix[localI][localI] += tetBuffer[vertexI];
            }

            tpr.gradNiGradNj(tetBuffer);

            for (edgeI = 0; edgeI < tetPointRef::nEdges; edgeI++)
            {
                curEdge = curTet.tetEdge(edgeI);

                localI = globalToLocalBuffer[curEdge.start()];
                localJ = globalToLocalBuffer[curEdge.end()];

                // Accumulate the off-diagonal; folding the matrix.
                // It is unknown if the localI or localJ is lower:
                // that depends on how the local points have been
                // selected.  Therefore, the matrix (which is
                // symmetric!) needs to be "folded"

                minIJ = min(localI, localJ);
                maxIJ = max(localI, localJ);

                // Take care of the directionality of the edge
                if (curEdge.start() < curEdge.end())
                {
                    denseMatrix[minIJ][maxIJ] += tetBuffer[edgeI];
                }
                else
                {
                    denseMatrix[minIJ][maxIJ] += tetBuffer[edgeI].T();
                }
            }
        }
    }
}


// Fill buffer with the volume integral distributed into vertices
void tetPolyMeshCellDecomp::volIntegral
(
    const label cellID,
    scalarField& buffer,
    const labelList& globalToLocalBuffer
) const
{
    // The addressing and volume distribution has been done tet-by-tet
    // and the combination is done later for the whole cell.

    const unallocFaceList& meshFaces = mesh_.faces();
    const unallocLabelList& owner = mesh_.faceOwner();
    const labelList& cellFaces = mesh_.cells()[cellID];

    const pointField& points = mesh_.points();
    const pointField& cellCentres = mesh_.cellCentres();

    label vertexI, localI;

    forAll (cellFaces, faceI)
    {
        const label curFaceID = cellFaces[faceI];

        face f = meshFaces[curFaceID];

        // Take care of owner/neighbour face reversal
        if (cellID == owner[curFaceID])
        {
            f = f.reverseFace();
        }

        // Grab the point zero of the face and create triangles by grabbing the
        // labels around the perimeter
        label pointZeroID = f[0];

        for (label i = 1; i < f.size() - 1; i++)
        {
            // Face is assumed to be inward-pointing
            tetCell curTet
            (
                pointZeroID,
                f[i],
                f[i + 1],
                cellOffset() + cellID
            );

            scalar quarterVolume =
                0.25*
                tetPointRef
                (
                    points[pointZeroID],
                    points[f[i]],
                    points[f[i + 1]],
                    cellCentres[cellID]
                ).mag();

            for
            (
                vertexI = 0;
                vertexI < tetPointRef::nVertices;
                vertexI++
            )
            {
                localI = globalToLocalBuffer[curTet[vertexI]];

                // Accumulate the diagonal
                buffer[localI] += quarterVolume;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
