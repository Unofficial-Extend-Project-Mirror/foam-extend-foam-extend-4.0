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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledPlane::createGeometry()
{
    clearGeom();
    faces_.clear();
    meshCells_.clear();

    if (triangulate())
    {
        // Count triangles
        label nTri = 0;
        forAll(cuttingPlane::faces(), faceI)
        {
            const face& f = cuttingPlane::faces()[faceI];

            nTri += f.nTriangles(points());
        }

        faces_.setSize(nTri);
        meshCells_.setSize(nTri);

        // split and fill mesh cell references
        nTri = 0;
        forAll(cuttingPlane::faces(), faceI)
        {
            const face& f = cuttingPlane::faces()[faceI];
            label cellId  = cuttingPlane::cells()[faceI];

            label fillIndex = nTri;

            f.triangles(points(), nTri, faces_);
            while (fillIndex < nTri)
            {
                meshCells_[fillIndex++] = cellId;
            }
        }
    }

    if (debug)
    {
        print(Pout);
        Pout << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledPlane::sampledPlane
(
    const word& name,
    const polyMesh& mesh,
    const plane& planeDesc,
    const word& zoneName,
    const bool triangulate
)
:
    sampledSurface(name, mesh, triangulate),
    cuttingPlane(planeDesc),
    zoneName_(zoneName),
    faces_(0),
    meshCells_(0)
{
    label zoneId = -1;
    if (zoneName_.size())
    {
        zoneId = mesh.cellZones().findZoneID(zoneName_);
        if (debug && zoneId < 0)
        {
            Info<< "cellZone \"" << zoneName_
                << "\" not found - using entire mesh"
                << endl;
        }
    }

    if (zoneId < 0)
    {
        reCut(mesh);
    }
    else
    {
        reCut(mesh, mesh.cellZones()[zoneId]);
    }

    createGeometry();
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
    faces_(0),
    meshCells_(0)
{

    // make plane relative to the coordinateSystem (Cartesian)
    if (dict.found("coordinateSystem"))
    {
        coordinateSystem cs(dict.subDict("coordinateSystem"));

        point  base = cs.globalPosition(planeDesc().refPoint());
        vector norm = cs.globalVector(planeDesc().normal());

        // assign the plane description
        static_cast<plane&>(*this) = plane(base, norm);
    }


    label zoneId = -1;
    if (dict.found("zone"))
    {
        dict.lookup("zone") >> zoneName_;
        zoneId = mesh.cellZones().findZoneID(zoneName_);
        if (debug && zoneId < 0)
        {
            Info<< "cellZone \"" << zoneName_
                << "\" not found - using entire mesh"
                << endl;
        }
    }


    if (zoneId < 0)
    {
        reCut(mesh);
    }
    else
    {
        reCut(mesh, mesh.cellZones()[zoneId]);
    }

    createGeometry();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledPlane::~sampledPlane()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sampledPlane::correct(const bool meshChanged)
{
    // Only change of mesh changes plane - zone restriction gets lost
    if (meshChanged)
    {
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

        createGeometry();
    }
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
