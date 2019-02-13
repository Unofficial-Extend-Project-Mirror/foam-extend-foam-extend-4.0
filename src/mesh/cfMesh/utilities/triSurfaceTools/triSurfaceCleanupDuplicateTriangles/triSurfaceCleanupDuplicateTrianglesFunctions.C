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

#include "triSurfaceCleanupDuplicateTriangles.H"
#include "triSurfModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceCleanupDuplicateTriangles::checkDuplicateTriangles()
{
    labelLongList newTriangleLabel(surf_.size(), -1);

    const VRWGraph& pointTriangles = surf_.pointFacets();

    //- check if there exist duplicate triangles
    label counter(0);

    forAll(surf_, triI)
    {
        if( newTriangleLabel[triI] != -1 )
            continue;

        newTriangleLabel[triI] = counter;
        ++counter;

        const labelledTri& tri = surf_[triI];

        forAll(pointTriangles[tri[0]], ptI)
        {
            const label triJ = pointTriangles(tri[0], ptI);

            if( triJ <= triI )
                continue;

            const labelledTri& otherTri = surf_[triJ];

            if( tri == otherTri )
                newTriangleLabel[triJ] = newTriangleLabel[triI];
        }
    }

    Info << "Found " << (newTriangleLabel.size()-counter)
        << " duplicate triangles" << endl;

    //- return if there exist no duplicate triangles
    if( counter == newTriangleLabel.size() )
        return;

    Info << "Current number of triangles" << surf_.size() << endl;
    Info << "New number of triangles " << counter << endl;

    //- create new list of triangles and store it in the surface mesh
    LongList<labelledTri> newTriangles(counter);

    forAll(newTriangleLabel, triI)
    {
        newTriangles[newTriangleLabel[triI]] = surf_[triI];
    }

    triSurfModifier(surf_).facetsAccess().transfer(newTriangles);
    surf_.updateFacetsSubsets(newTriangleLabel);

    surf_.clearAddressing();
    surf_.clearGeometry();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
