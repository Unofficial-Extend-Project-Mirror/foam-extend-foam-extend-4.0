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

\*---------------------------------------------------------------------------*/

#include "polyMeshGenGeometryModification.H"
#include "dictionary.H"

namespace Foam
{

// * * * * * * * * * * * * * * Private member functions* * * * * * * * * * * //

void polyMeshGenGeometryModification::checkModification()
{
    if( meshDict_.found("anisotropicSources") )
    {
        modificationActive_ = true;

        const dictionary& anisotropicDict =
            meshDict_.subDict("anisotropicSources");

        coordinateModifierPtr_ = new coordinateModifier(anisotropicDict);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

polyMeshGenGeometryModification::polyMeshGenGeometryModification
(
    polyMeshGen& mesh,
    const dictionary& meshDict
)
:
    mesh_(mesh),
    meshDict_(meshDict),
    coordinateModifierPtr_(nullptr),
    modificationActive_(false)
{
    checkModification();
}

polyMeshGenGeometryModification::~polyMeshGenGeometryModification()
{
    deleteDemandDrivenData(coordinateModifierPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool polyMeshGenGeometryModification::activeModification() const
{
    return modificationActive_;
}

void polyMeshGenGeometryModification::modifyGeometry()
{
    if( !modificationActive_ )
    {
        WarningIn
        (
            "const triSurf* polyMeshGenGeometryModification"
            "::modifyGeometry() const"
        ) << "Modification is not active" << endl;

        return;
    }

    pointFieldPMG& pts = mesh_.points();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(pts, pointI)
        pts[pointI] = coordinateModifierPtr_->modifiedPoint(pts[pointI]);
}

void polyMeshGenGeometryModification::revertGeometryModification()
{
    if( !modificationActive_ )
    {
        WarningIn
        (
            "const triSurf* polyMeshGenGeometryModification"
            "::revertGeometryModification() const"
        ) << "Modification is not active" << endl;

        return;
    }

    pointFieldPMG& pts = mesh_.points();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(pts, pointI)
        pts[pointI] =
            coordinateModifierPtr_->backwardModifiedPoint(pts[pointI]);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
