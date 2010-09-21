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

#include "tetPolyMeshFaceDecomp.H"
#include "tetPointRef.H"
#include "tetPointRef.H"
#include "SquareMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label tetPolyMeshFaceDecomp::nPoints() const
{
    if (nPoints_ < 0)
    {
        nPoints_ =
            mesh_.nPoints()
          + mesh_.nFaces()
          + mesh_.nCells();
    }

    return nPoints_;
}


label tetPolyMeshFaceDecomp::nEdges() const
{
    if (nEdges_ < 0)
    {
        nEdges_ = mesh_.nEdges();

        const labelListList& pf = mesh_.pointFaces();

        const labelListList& pc = mesh_.pointCells();

        forAll (pf, pointI)
        {
            nEdges_ += pf[pointI].size();

            nEdges_ += pc[pointI].size();
        }

        const cellList& cellFaces = mesh_.cells();

        forAll (cellFaces, cellI)
        {
            nEdges_ += cellFaces[cellI].size();
        }
    }

    return nEdges_;
}


label tetPolyMeshFaceDecomp::nTets() const
{
    if (nTets_ < 0)
    {
        // count the cells
        nTets_ = 0;

        const cellList& polyCells = mesh_.cells();

        forAll (polyCells, cellI)
        {
            nTets_ += nTetsForCell(cellI);
        }
    }

    return nTets_;
}


// Return number of tetrahedra in decomposition for cell
label tetPolyMeshFaceDecomp::nTetsForCell(const label cellID) const
{
    const unallocFaceList& f = mesh_.faces();

    label nTetrasForCell = 0;

    const labelList& cellFaces = mesh_.cells()[cellID];

    forAll (cellFaces, faceI)
    {
        nTetrasForCell += f[cellFaces[faceI]].nEdges();
    }

    return nTetrasForCell;
}


tmp<pointField> tetPolyMeshFaceDecomp::points() const
{
    tmp<pointField> ttetPoints(new pointField(nPoints()));
    pointField& tetPoints = ttetPoints();

    const pointField& points = mesh_.points();

    const pointField& faceCentres = mesh_.faceCentres();

    const pointField& cellCentres = mesh_.cellCentres();

    label tetPointI = 0;

    forAll (points, pointI)
    {
        tetPoints[tetPointI] = points[pointI];
        tetPointI++;
    }

    forAll (faceCentres, faceI)
    {
        tetPoints[tetPointI] = faceCentres[faceI];
        tetPointI++;
    }

    forAll (cellCentres, cellI)
    {
        tetPoints[tetPointI] = cellCentres[cellI];
        tetPointI++;
    }

    return ttetPoints;
}


cellShapeList tetPolyMeshFaceDecomp::tetCells() const
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
tetCellList tetPolyMeshFaceDecomp::tets(const label cellID) const
{
    const unallocFaceList& f = mesh_.faces();

    const unallocLabelList& owner = mesh_.faceOwner();

    // Initialise the size of the return
    tetCellList t(nTetsForCell(cellID));

    label nTetras = 0;

    const labelList& cellFaces = mesh_.cells()[cellID];

    // A tet is created from the face edge, face centre and cell centre
    forAll (cellFaces, faceI)
    {
        const label curFace = cellFaces[faceI];

        edgeList faceEdges;

        // Take care of owner/neighbour face reversal
        if (cellID == owner[curFace])
        {
            faceEdges = f[curFace].reverseFace().edges();
        }
        else
        {
            faceEdges = f[curFace].edges();
        }

        forAll (faceEdges, edgeI)
        {
            // Face is assumed to be inward-pointing
            t[nTetras] =
                tetCell
                (
                    faceEdges[edgeI].start(),
                    faceEdges[edgeI].end(),
                    curFace + faceOffset_,
                    cellID + cellOffset_
                );

            nTetras++;
        }
    }

    return t;
}


void tetPolyMeshFaceDecomp::gradNiDotGradNj
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
    const pointField& faceCentres = mesh_.faceCentres();
    const pointField& cellCentres = mesh_.cellCentres();

    scalarField tetBuffer(tetPointRef::nEdges);

    // For optimisation, point allocation is taken out of the loop
    label vertexI, edgeI, localI, localJ;
    edge curEdge;

    forAll (cellFaces, faceI)
    {
        const label curFaceID = cellFaces[faceI];

        edgeList faceEdges;

        if (cellID == owner[curFaceID])
        {
            faceEdges = meshFaces[curFaceID].reverseFace().edges();
        }
        else
        {
            faceEdges = meshFaces[curFaceID].edges();
        }

        forAll (faceEdges, i)
        {
            tetCell curTet
            (
                faceEdges[i].start(),
                faceEdges[i].end(),
                faceOffset() + curFaceID,
                cellOffset() + cellID
            );

            tetPointRef tpr
            (
                points[faceEdges[i].start()],
                points[faceEdges[i].end()],
                faceCentres[curFaceID],
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

            tpr.gradNiDotGradNj(tetBuffer);

            for (edgeI = 0; edgeI < tetPointRef::nEdges; edgeI++)
            {
                curEdge = curTet.tetEdge(edgeI);

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


void tetPolyMeshFaceDecomp::gradNiGradNj
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
    const pointField& faceCentres = mesh_.faceCentres();
    const pointField& cellCentres = mesh_.cellCentres();

    tensorField tetBuffer(tetPointRef::nEdges);

    // For optimisation, point allocation is taken out of the loop
    label vertexI, edgeI, localI, localJ;
    edge curEdge;

    forAll (cellFaces, faceI)
    {
        const label curFaceID = cellFaces[faceI];

        edgeList faceEdges;

        if (cellID == owner[curFaceID])
        {
            faceEdges = meshFaces[curFaceID].reverseFace().edges();
        }
        else
        {
            faceEdges = meshFaces[curFaceID].edges();
        }

        forAll (faceEdges, i)
        {
            tetCell curTet
            (
                faceEdges[i].start(),
                faceEdges[i].end(),
                faceOffset() + curFaceID,
                cellOffset() + cellID
            );

            tetPointRef tpr
            (
                points[faceEdges[i].start()],
                points[faceEdges[i].end()],
                faceCentres[curFaceID],
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

                // Add coefficient and its transpose to upper and lower
                // triangle.  HJ, 4/Dec/2006
                denseMatrix[localI][localJ] += tetBuffer[edgeI];
                denseMatrix[localJ][localI] += tetBuffer[edgeI].T();
            }
        }
    }
}


// Fill buffer with the volume integral distributed into vertices
void tetPolyMeshFaceDecomp::volIntegral
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
    const pointField& faceCentres = mesh_.faceCentres();
    const pointField& cellCentres = mesh_.cellCentres();

    label vertexI, localI;

    forAll (cellFaces, faceI)
    {
        const label curFaceID = cellFaces[faceI];

        edgeList faceEdges;

        // Take care of owner/neighbour face reversal
        if (cellID == owner[curFaceID])
        {
            faceEdges = meshFaces[curFaceID].reverseFace().edges();
        }
        else
        {
            faceEdges = meshFaces[curFaceID].edges();
        }

        forAll (faceEdges, edgeI)
        {
            // Face is assumed to be inward-pointing
            tetCell curTet
            (
                faceEdges[edgeI].start(),
                faceEdges[edgeI].end(),
                faceOffset() + curFaceID,
                cellOffset() + cellID
            );

            scalar quarterVolume =
                0.25*
                tetPointRef
                (
                    points[faceEdges[edgeI].start()],
                    points[faceEdges[edgeI].end()],
                    faceCentres[curFaceID],
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
