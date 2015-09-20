/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

InClass
    vtkPV3Foam

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3FoamPoints_H
#define vtkPV3FoamPoints_H

// VTK includes
#include "vtkPoints.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline void vtkPV3FoamInsertNextPoint
(
    vtkPoints *points,
    const Foam::point& p
)
{
    points->InsertNextPoint(p.x(), p.y(), p.z());
}

#if 0
// this should be faster, but didn't get it working ...
inline void vtkPV3FoamSetPoint
(
    vtkPoints *points,
    const Foam::label id,
    const Foam::point& p
)
{
    points->SetPoint(id, p.x(), p.y(), p.z());
}


// Convert Foam mesh vertices to VTK
inline vtkPoints* vtkPV3FoamVTKPoints(const Foam::pointField& points)
{
    vtkPoints *vtkpoints = vtkPoints::New();
    vtkpoints->SetNumberOfPoints(points.size());
    forAll(points, i)
    {
        const Foam::point& p = points[i];
        vtkpoints->SetPoint(i, p.x(), p.y(), p.z());
    }

    return vtkpoints;
}

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
