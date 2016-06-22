/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

\*---------------------------------------------------------------------------*/

#include "MeshedSurfaceIOAllocator.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MeshedSurfaceIOAllocator::MeshedSurfaceIOAllocator
(
    const IOobject& ioPoints,
    const IOobject& ioFaces,
    const IOobject& ioZones
)
:
    points_(ioPoints),
    faces_(ioFaces),
    zones_(ioZones)
{}


Foam::MeshedSurfaceIOAllocator::MeshedSurfaceIOAllocator
(
    const IOobject& ioPoints,
    const pointField& points,
    const IOobject& ioFaces,
    const faceList& faces,
    const IOobject& ioZones,
    const surfZoneList& zones
)
:
    points_(ioPoints, points),
    faces_(ioFaces, faces),
    zones_(ioZones, zones)
{}


Foam::MeshedSurfaceIOAllocator::MeshedSurfaceIOAllocator
(
    const IOobject& ioPoints,
    const Xfer< pointField >& points,
    const IOobject& ioFaces,
    const Xfer< faceList >& faces,
    const IOobject& ioZones,
    const Xfer< surfZoneList >& zones
)
:
    points_(ioPoints, points),
    faces_(ioFaces, faces),
    zones_(ioZones, zones)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MeshedSurfaceIOAllocator::clear()
{
    points_.clear();
    faces_.clear();
    zones_.clear();
}


void Foam::MeshedSurfaceIOAllocator::resetFaces
(
    const Xfer<List<face> >& faces,
    const Xfer<surfZoneList>& zones
)
{
    if (!faces().empty())
    {
        faces_.transfer(faces());
    }

    if (!zones().empty())
    {
        zones_.transfer(zones());
    }
}


void Foam::MeshedSurfaceIOAllocator::reset
(
    const Xfer< pointField >& points,
    const Xfer< faceList >& faces,
    const Xfer< surfZoneList >& zones
)
{
    // Take over new primitive data.
    // Optimized to avoid overwriting data at all
    if (!points().empty())
    {
        points_.transfer(points());
    }

    resetFaces(faces, zones);
}


void Foam::MeshedSurfaceIOAllocator::reset
(
    const Xfer< List<point> >& points,
    const Xfer< faceList >& faces,
    const Xfer< surfZoneList >& zones
)
{
    // Take over new primitive data.
    // Optimized to avoid overwriting data at all
    if (!points().empty())
    {
        points_.transfer(points());
    }

    resetFaces(faces, zones);
}


// ************************************************************************* //
