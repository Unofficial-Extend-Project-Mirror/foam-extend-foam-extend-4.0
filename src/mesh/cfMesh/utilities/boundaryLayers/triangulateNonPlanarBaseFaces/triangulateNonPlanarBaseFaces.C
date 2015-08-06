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

#include "triangulateNonPlanarBaseFaces.H"
#include "dictionary.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

triangulateNonPlanarBaseFaces::triangulateNonPlanarBaseFaces
(
    polyMeshGen& mesh
)
:
    mesh_(mesh),
    invertedCell_(mesh_.cells().size(), false),
    decomposeFace_(mesh_.faces().size(), false),
    tol_(0.5)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

triangulateNonPlanarBaseFaces::~triangulateNonPlanarBaseFaces()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triangulateNonPlanarBaseFaces::setRelativeTolerance(const scalar tol)
{
    tol_ = tol;
}

void triangulateNonPlanarBaseFaces::triangulateLayers()
{
    if( findNonPlanarBoundaryFaces() )
    {
        Info << "Decomposing twisted boundary faces" << endl;

        decomposeBoundaryFaces();

        decomposeCellsIntoPyramids();
    }
    else
    {
        Info << "All boundary faces are flat" << endl;
    }
}

void triangulateNonPlanarBaseFaces::readSettings
(
    const dictionary& meshDict,
    triangulateNonPlanarBaseFaces& triangulator
)
{
    if( meshDict.found("boundaryLayers") )
    {
        const dictionary& layersDict = meshDict.subDict("boundaryLayers");

        if( layersDict.found("optimisationParameters") )
        {
            const dictionary& optLayerDict =
                layersDict.subDict("optimisationParameters");

            if( optLayerDict.found("relFlatnessTol") )
            {
                const scalar relTol =
                    readScalar(optLayerDict.lookup("relFlatnessTol"));

                triangulator.setRelativeTolerance(relTol);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
