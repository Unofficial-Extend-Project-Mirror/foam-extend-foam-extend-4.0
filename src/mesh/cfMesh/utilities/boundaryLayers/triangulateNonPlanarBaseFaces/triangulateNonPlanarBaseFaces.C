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
