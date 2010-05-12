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

#include "faceDecompIsoSurfaceCuts.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceDecompIsoSurfaceCuts::constructEdgeCuts
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
    labelHashSet cutFaceCentres(nBFaces);
    labelHashSet cutCellCentres(nBFaces);

    DynamicList<label> meshCutEdges(nBFaces);
    DynamicList<scalar> meshEdgeWeights(nBFaces);

    DynamicList<pyramidEdge> cutPyrEdges(nBFaces);
    DynamicList<scalar> pyrEdgeWeights(nBFaces);

    DynamicList<faceEdge> cutFaceEdges(nBFaces);
    DynamicList<scalar> faceEdgeWeights(nBFaces);

    DynamicList<centreEdge> cutCentreEdges(nBFaces);
    DynamicList<scalar> centreEdgeWeights(nBFaces);


    // Walk over faces
    forAll(mesh.faces(), faceI)
    {
        const face& f = mesh.faces()[faceI];

        // Face centre value
        scalar fcVal = 0;

        forAll(f, fp)
        {
            fcVal += ptField[f[fp]];
        }
        fcVal /= f.size();


        label ownerI = mesh.faceOwner()[faceI];
        label neighbourI = -1;
        if (mesh.isInternalFace(faceI))
        {
            neighbourI = mesh.faceNeighbour()[faceI];
        }

        scalar weight;

        //
        // Check faceCentre-cellCentre crossing for owner and neigbour
        //

        if (crosses(isoVal, fcVal, volField[ownerI], weight))
        {
            mark(ownerI, cutCells);

            if (weight < tol)
            {
                mark(faceI, cutFaceCentres);
            }
            else if (weight > 1-tol)
            {
                mark(ownerI, cutCellCentres);
            }
            else
            {
                cutCentreEdges.append(centreEdge(faceI, true));
                centreEdgeWeights.append(weight);
            }
        }

        if (neighbourI != -1)
        {
            if (crosses(isoVal, fcVal, volField[neighbourI], weight))
            {
                mark(neighbourI, cutCells);

                if (weight < tol)
                {
                    mark(faceI, cutFaceCentres);
                }
                else if (weight > 1-tol)
                {
                    mark(neighbourI, cutCellCentres);
                }
                else
                {
                    cutCentreEdges.append(centreEdge(faceI, false));
                    centreEdgeWeights.append(weight);
                }
            }
        }


        forAll(f, fp)
        {
            //
            // Check face internal (faceEdge) edge crossing
            //

            if (crosses(isoVal, ptField[f[fp]], fcVal, weight))
            {
                mark(ownerI, cutCells);
                if (neighbourI != -1)
                {
                    mark(neighbourI, cutCells);
                }

                if (weight < tol)
                {
                    mark(f[fp], cutVerts);
                }
                else if (weight > 1-tol)
                {
                    mark(faceI, cutFaceCentres);
                }
                else
                {
                    cutFaceEdges.append(faceEdge(faceI, fp));
                    faceEdgeWeights.append(weight);
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

    cutFaceEdges.shrink();
    faceEdgeWeights.shrink();

    cutCentreEdges.shrink();
    centreEdgeWeights.shrink();



    // Tranfer lists to faceDecompCuts

    cells_ = cutCells.toc();
    meshVerts_ = cutVerts.toc();
    meshFaceCentres_ = cutFaceCentres.toc();
    meshCellCentres_ = cutCellCentres.toc();

    meshEdges_.transfer(meshCutEdges);
    meshEdgeWeights_.transfer(meshEdgeWeights);

    pyrEdges_.transfer(cutPyrEdges);
    pyrEdgeWeights_.transfer(pyrEdgeWeights);

    faceEdges_.transfer(cutFaceEdges);
    faceEdgeWeights_.transfer(faceEdgeWeights);

    centreEdges_.transfer(cutCentreEdges);
    centreEdgeWeights_.transfer(centreEdgeWeights);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceDecompIsoSurfaceCuts::faceDecompIsoSurfaceCuts
(
    const volScalarField& volField,
    const pointScalarField& ptField,
    const scalar isoVal,
    const scalar tol
)
:
    faceDecompCuts
    (
        volField.mesh(),
        labelList(),

        labelList(),
        labelList(),
        labelList(),

        labelList(),            // mesh edges
        scalarField(),

        List<pyramidEdge>(),    // pyramid edges
        scalarField(),

        List<faceEdge>(),       // face edges
        scalarField(),

        List<centreEdge>(),     // cell centre edges
        scalarField()
    )
{
    constructEdgeCuts(volField, ptField, isoVal, tol);
}


// Construct from interpolator and field (does interpolation)
Foam::faceDecompIsoSurfaceCuts::faceDecompIsoSurfaceCuts
(
    const volScalarField& volField,
    const volPointInterpolation& pInterp,
    const scalar isoVal,
    const scalar tol
)
:
    faceDecompCuts
    (
        volField.mesh(),
        labelList(),

        labelList(),
        labelList(),
        labelList(),

        labelList(),            // mesh edges
        scalarField(),

        List<pyramidEdge>(),    // pyramid edges
        scalarField(),

        List<faceEdge>(),       // face edges
        scalarField(),

        List<centreEdge>(),     // cell centre edges
        scalarField()
    )
{
    // Get field on vertices
    tmp<pointScalarField> tptField = pInterp.interpolate(volField);
    const pointScalarField& ptField = tptField();


    constructEdgeCuts(volField, ptField, isoVal, tol);
}


// ************************************************************************* //
