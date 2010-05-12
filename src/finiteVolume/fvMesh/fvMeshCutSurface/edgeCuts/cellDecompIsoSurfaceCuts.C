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

Description

\*---------------------------------------------------------------------------*/

#include "cellDecompIsoSurfaceCuts.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellDecompIsoSurfaceCuts::constructEdgeCuts
(
    const volScalarField& volField,
    const pointScalarField& ptField,
    const scalar isoVal,
    const scalar tol
)
{
    const primitiveMesh& mesh = volField.mesh();

    // Intermediate storage for labels of cut edges and cut points
    label nBFaces = mesh.nFaces() - mesh.nInternalFaces();

    labelHashSet cutCells(nBFaces);
    labelHashSet cutVerts(nBFaces);
    labelHashSet cutCellCentres(nBFaces);

    DynamicList<label> meshCutEdges(nBFaces);
    DynamicList<scalar> meshEdgeWeights(nBFaces);

    DynamicList<pyramidEdge> cutPyrEdges(nBFaces);
    DynamicList<scalar> pyrEdgeWeights(nBFaces);

    DynamicList<diagonalEdge> cutDiagEdges(nBFaces);
    DynamicList<scalar> diagEdgeWeights(nBFaces);


    // Walk over faces
    forAll(mesh.faces(), faceI)
    {
        const face& f = mesh.faces()[faceI];

        label ownerI = mesh.faceOwner()[faceI];
        label neighbourI = -1;
        if (mesh.isInternalFace(faceI))
        {
            neighbourI = mesh.faceNeighbour()[faceI];
        }

        scalar weight;

        //
        // Check face diagonal. Simple triangulation from f[0]
        //

        for(label fp = 2; fp < f.size() - 1; fp++)
        {
            if (crosses(isoVal, ptField[f[0]], ptField[f[fp]], weight))
            {
                mark(ownerI, cutCells);
                if (neighbourI != -1)
                {
                    mark(neighbourI, cutCells);
                }

                if (weight < tol)
                {
                    mark(f[0], cutVerts);
                }
                else if (weight > 1-tol)
                {
                    mark(f[fp], cutVerts);
                }
                else
                {
                    cutDiagEdges.append(diagonalEdge(faceI, 0, fp));
                    diagEdgeWeights.append(weight);
                }
            }
        }
    }
    

    // Walk over all mesh points
    forAll(mesh.points(), pointI)
    {
        const labelList& myCells = mesh.pointCells()[pointI];

        forAll(myCells, myCellI)
        {
            label cellI = myCells[myCellI];

            //
            // Check pyramid edge crossing
            //

            scalar weight;

            if (crosses(isoVal, ptField[pointI], volField[cellI], weight))
            {
                mark(cellI, cutCells);

                if (weight < tol)
                {
                    mark(pointI, cutVerts);
                }
                else if (weight > 1-tol)
                {
                    mark(cellI, cutCellCentres);
                }
                else
                {
                    cutPyrEdges.append(pyramidEdge(pointI, cellI));
                    pyrEdgeWeights.append(weight);
                }
            }
        }
    }            


    // Walk over all mesh edges
    forAll(mesh.edges(), edgeI)
    {
        const edge& e = mesh.edges()[edgeI];

        scalar weight;

        if (crosses(isoVal, ptField[e.start()], ptField[e.end()], weight))
        {
            mark(mesh.edgeCells()[edgeI], cutCells);

            if (weight < tol)
            {
                mark(e.start(), cutVerts);
            }
            else if (weight > 1-tol)
            {
                mark(e.end(), cutVerts);
            }
            else
            {
                meshCutEdges.append(edgeI);
                meshEdgeWeights.append(weight);
            }
        }
    }

    meshCutEdges.shrink();
    meshEdgeWeights.shrink();

    cutPyrEdges.shrink();
    pyrEdgeWeights.shrink();

    cutDiagEdges.shrink();
    diagEdgeWeights.shrink();


    // Tranfer lists to cellDecompCuts

    cells_ = cutCells.toc();

    meshVerts_ = cutVerts.toc();
    meshCellCentres_ = cutCellCentres.toc();

    meshEdges_.transfer(meshCutEdges);
    meshEdgeWeights_.transfer(meshEdgeWeights);

    pyrEdges_.transfer(cutPyrEdges);
    pyrEdgeWeights_.transfer(pyrEdgeWeights);

    diagEdges_.transfer(cutDiagEdges);
    diagEdgeWeights_.transfer(diagEdgeWeights);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellDecompIsoSurfaceCuts::cellDecompIsoSurfaceCuts
(
    const volScalarField& volField,
    const pointScalarField& ptField,
    const scalar isoVal,
    const scalar tol
)
:
    cellDecompCuts
    (
        volField.mesh(),
        labelList(),

        labelList(),            // mesh vertices
        labelList(),            // mesh cell centres

        labelList(),            // mesh edges
        scalarField(),

        List<pyramidEdge>(),    // pyramid edges
        scalarField(),

        List<diagonalEdge>(),   // face diagonal edges
        scalarField()
    )
{
    constructEdgeCuts(volField, ptField, isoVal, tol);
}


// Construct from interpolator and field (does interpolation)
Foam::cellDecompIsoSurfaceCuts::cellDecompIsoSurfaceCuts
(
    const volScalarField& volField,
    const volPointInterpolation& pInterp,
    const scalar isoVal,
    const scalar tol
)
:
    cellDecompCuts
    (
        volField.mesh(),
        labelList(),

        labelList(),            // mesh vertices
        labelList(),            // mesh cell centres

        labelList(),            // mesh edges
        scalarField(),

        List<pyramidEdge>(),    // pyramid edges
        scalarField(),

        List<diagonalEdge>(),   // face diagonal edges
        scalarField()
    )
{
    // Get field on vertices
    tmp<pointScalarField> tptField = pInterp.interpolate(volField);
    const pointScalarField& ptField = tptField();


    constructEdgeCuts(volField, ptField, isoVal, tol);
}


// ************************************************************************* //
