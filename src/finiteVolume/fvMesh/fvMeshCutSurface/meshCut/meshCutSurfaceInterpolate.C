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

#include "meshCutSurface.H"
#include "pyramidEdge.H"
#include "faceEdge.H"
#include "centreEdge.H"
#include "diagonalEdge.H"
#include "faceDecompCuts.H"
#include "cellDecompCuts.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class T>
Foam::tmp<Foam::Field<T> > Foam::meshCutSurface::interpolate
(
    const faceDecompCuts& cuts,
    const Field<T>& vField,
    const Field<T>& fField,
    const Field<T>& pField
)
{
    const primitiveMesh& mesh = cuts.mesh();

    tmp<Field<T> > tvals(new Field<T>(cuts.size()));
    Field<T>& vals = tvals();

    label pointI = 0;

    forAll(cuts.meshVerts(), cutVertI)
    {
        vals[pointI++] = pField[cuts.meshVerts()[cutVertI]];
    }

    forAll(cuts.meshFaceCentres(), cutFaceCentreI)
    {
        label faceI = cuts.meshFaceCentres()[cutFaceCentreI];

        vals[pointI++] = fField[faceI];
    }

    forAll(cuts.meshCellCentres(), cutCellCentreI)
    {
        label cellI = cuts.meshCellCentres()[cutCellCentreI];

        vals[pointI++] = vField[cellI];
    }

    forAll(cuts.meshEdges(), meshCutSurfaceEdgeI)
    {
        label edgeI = cuts.meshEdges()[meshCutSurfaceEdgeI];

        const edge& e = mesh.edges()[edgeI];

        scalar weight = cuts.meshEdgeWeights()[meshCutSurfaceEdgeI];

        vals[pointI++] = weight*pField[e.start()] + (1-weight)*pField[e.end()];
    }

    forAll(cuts.pyrEdges(), edgeI)
    {
        const pyramidEdge& e = cuts.pyrEdges()[edgeI];

        scalar weight = cuts.pyrEdgeWeights()[edgeI];

        vals[pointI++] = e.interpolate(vField, pField, weight);
    }

    forAll(cuts.centreEdges(), edgeI)
    {
        const centreEdge& e = cuts.centreEdges()[edgeI];

        scalar weight = cuts.centreEdgeWeights()[edgeI];

        vals[pointI++] = e.interpolate(mesh, vField, fField, weight);
    }

    forAll(cuts.faceEdges(), edgeI)
    {
        const faceEdge& e = cuts.faceEdges()[edgeI];

        scalar weight = cuts.faceEdgeWeights()[edgeI];

        vals[pointI++] = e.interpolate(mesh, fField, pField, weight);
    }

    return tvals;
}


template <class T>
Foam::tmp<Foam::Field<T> > Foam::meshCutSurface::interpolate
(
    const cellDecompCuts& cuts,
    const Field<T>& vField,
    const Field<T>& pField
)
{
    const primitiveMesh& mesh = cuts.mesh();

    tmp<Field<T> > tvals(new Field<T>(cuts.size()));
    Field<T>& vals = tvals();

    label pointI = 0;

    forAll(cuts.meshVerts(), cutVertI)
    {
        vals[pointI++] = pField[cuts.meshVerts()[cutVertI]];
    }

    forAll(cuts.meshCellCentres(), cutCellCentreI)
    {
        label cellI = cuts.meshCellCentres()[cutCellCentreI];

        vals[pointI++] = vField[cellI];
    }

    forAll(cuts.meshEdges(), meshCutSurfaceEdgeI)
    {
        label edgeI = cuts.meshEdges()[meshCutSurfaceEdgeI];

        const edge& e = mesh.edges()[edgeI];

        scalar weight = cuts.meshEdgeWeights()[meshCutSurfaceEdgeI];

        vals[pointI++] = weight*pField[e.start()] + (1-weight)*pField[e.end()];
    }

    forAll(cuts.pyrEdges(), edgeI)
    {
        const pyramidEdge& e = cuts.pyrEdges()[edgeI];

        scalar weight = cuts.pyrEdgeWeights()[edgeI];

        vals[pointI++] = e.interpolate(vField, pField, weight);
    }

    forAll(cuts.diagEdges(), edgeI)
    {
        const diagonalEdge& e = cuts.diagEdges()[edgeI];

        scalar weight = cuts.diagEdgeWeights()[edgeI];

        vals[pointI++] = e.interpolate(mesh, pField, weight);
    }

    return tvals;
}


// ************************************************************************* //
