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

\*---------------------------------------------------------------------------*/

#include "sampledPlane.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "volFields.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledPlane, 0);
    addNamedToRunTimeSelectionTable(sampledSurface, sampledPlane, word, plane);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledPlane::sampledPlane
(
    const word& name,
    const polyMesh& mesh,
    const plane& planeDesc,
    const word& zoneName
)
:
    sampledSurface(name, mesh),
    cuttingPlane(planeDesc),
    zoneName_(zoneName),
    needsUpdate_(true)
{
    if (debug && zoneName_.size())
    {
        if (mesh.cellZones().findZoneID(zoneName_) < 0)
        {
            Info<< "cellZone \"" << zoneName_
                << "\" not found - using entire mesh" << endl;
        }
    }
}


Foam::sampledPlane::sampledPlane
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    cuttingPlane(plane(dict.lookup("basePoint"), dict.lookup("normalVector"))),
    zoneName_(word::null),
    needsUpdate_(true)
{
    // make plane relative to the coordinateSystem (Cartesian)
    // allow lookup from global coordinate systems
    if (dict.found("coordinateSystem"))
    {
        coordinateSystem cs(dict, mesh);

        point  base = cs.globalPosition(planeDesc().refPoint());
        vector norm = cs.globalVector(planeDesc().normal());

        // assign the plane description
        static_cast<plane&>(*this) = plane(base, norm);
    }

    dict.readIfPresent("zone", zoneName_);

    if (debug && zoneName_.size())
    {
        if (mesh.cellZones().findZoneID(zoneName_) < 0)
        {
            Info<< "cellZone \"" << zoneName_
                << "\" not found - using entire mesh" << endl;
        }
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledPlane::~sampledPlane()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledPlane::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledPlane::expire()
{
    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledPlane::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();

    label zoneId = -1;
    if (zoneName_.size())
    {
        zoneId = mesh().cellZones().findZoneID(zoneName_);
    }

    if (zoneId < 0)
    {
        reCut(mesh());
    }
    else
    {
        reCut(mesh(), mesh().cellZones()[zoneId]);
    }

    if (debug)
    {
        print(Pout);
        Pout << endl;
    }

    needsUpdate_ = false;
    return true;
}


Foam::tmp<Foam::scalarField>
Foam::sampledPlane::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField>
Foam::sampledPlane::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledPlane::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledPlane::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField>
Foam::sampledPlane::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField>
Foam::sampledPlane::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledPlane::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledPlane::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledPlane::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledPlane::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledPlane::print(Ostream& os) const
{
    os  << "sampledPlane: " << name() << " :"
        << "  base:" << refPoint()
        << "  normal:" << normal()
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
