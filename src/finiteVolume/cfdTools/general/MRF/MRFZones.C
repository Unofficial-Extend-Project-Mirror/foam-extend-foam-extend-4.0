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

#include "MRFZones.H"
#include "Time.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<MRFZone>, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFZones::MRFZones(const fvMesh& mesh)
:
    IOPtrList<MRFZone>
    (
        IOobject
        (
            "MRFZones",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        MRFZone::iNew(mesh)
    ),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::MRFZones::fluxCorrection() const
{
    tmp<surfaceScalarField> tMRFZonesPhiCorr
    (
        new surfaceScalarField
        (
            IOobject
            (
                "MRFZonesPhiCorr",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimVelocity*dimArea, 0)
        )
    );
    surfaceScalarField& MRFZonesPhiCorr = tMRFZonesPhiCorr();

    forAll(*this, i)
    {
        operator[](i).relativeFlux(MRFZonesPhiCorr);
    }

    return tMRFZonesPhiCorr;
}

Foam::tmp<Foam::surfaceVectorField> Foam::MRFZones::faceU() const
{
    tmp<surfaceVectorField> tMRFZonesFaceU
    (
        new surfaceVectorField
        (
            IOobject
            (
                "MRFZonesFaceU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimVelocity, vector::zero)
        )
    );
    surfaceVectorField& MRFZonesFaceU = tMRFZonesFaceU();

    forAll(*this, i)
    {
        operator[](i).faceU(MRFZonesFaceU);
    }

    return tMRFZonesFaceU;
}

void Foam::MRFZones::addCoriolis(fvVectorMatrix& UEqn) const
{
    forAll(*this, i)
    {
        operator[](i).addCoriolis(UEqn);
    }
}


void Foam::MRFZones::relativeFlux(surfaceScalarField& phi) const
{
    forAll(*this, i)
    {
        operator[](i).relativeFlux(phi);
    }
}


void Foam::MRFZones::absoluteFlux(surfaceScalarField& phi) const
{
    forAll(*this, i)
    {
        operator[](i).absoluteFlux(phi);
    }
}


void Foam::MRFZones::addCoriolis
(
    const volScalarField& rho,
    fvVectorMatrix& UEqn
) const
{
    forAll(*this, i)
    {
        operator[](i).addCoriolis(rho, UEqn);
    }
}


void Foam::MRFZones::relativeFlux
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    forAll(*this, i)
    {
        operator[](i).relativeFlux(rho, phi);
    }
}


void Foam::MRFZones::absoluteFlux
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    forAll(*this, i)
    {
        operator[](i).absoluteFlux(rho, phi);
    }
}


void Foam::MRFZones::relativeVelocity(volVectorField& U) const
{
    forAll(*this, i)
    {
        operator[](i).relativeVelocity(U);
    }
}


void Foam::MRFZones::absoluteVelocity(volVectorField& U) const
{
    forAll(*this, i)
    {
        operator[](i).absoluteVelocity(U);
    }
}


void Foam::MRFZones::correctBoundaryVelocity(volVectorField& U) const
{
    forAll(*this, i)
    {
        operator[](i).correctBoundaryVelocity(U);
    }
}


// ************************************************************************* //
