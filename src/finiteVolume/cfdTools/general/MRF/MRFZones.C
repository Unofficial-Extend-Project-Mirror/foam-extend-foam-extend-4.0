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
#include "fvc.H"

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

Foam::tmp<Foam::surfaceVectorField> Foam::MRFZones::meshPhi() const
{
    tmp<surfaceVectorField> tMRFZonesPhiCorr
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
    surfaceVectorField& MRFZonesPhiCorr = tMRFZonesPhiCorr();

    forAll(*this, i)
    {
        operator[](i).meshPhi(MRFZonesPhiCorr);
    }

    return tMRFZonesPhiCorr;
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


Foam::tmp<Foam::volScalarField> Foam::MRFZones::Su
(
    const volScalarField& phi
) const
{
    tmp<volScalarField> tPhiSource
    (
        new volScalarField
        (
            IOobject
            (
                phi.name() + "Source",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", phi.dimensions()/dimTime, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& source = tPhiSource();

    volVectorField gradPhi = fvc::grad(phi);

    forAll(*this, i)
    {
        operator[](i).Su(phi, gradPhi, source);
    }

    return tPhiSource;
}


Foam::tmp<Foam::volVectorField> Foam::MRFZones::Su
(
    const volVectorField& phi
) const
{
    tmp<volVectorField> tPhiSource
    (
        new volVectorField
        (
            IOobject
            (
                phi.name() + "Source",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", phi.dimensions()/dimTime, vector::zero),
            zeroGradientFvPatchVectorField::typeName
        )
    );
    volVectorField& source = tPhiSource();

    volTensorField gradPhi = fvc::grad(phi);

    forAll(*this, i)
    {
        operator[](i).Su(phi, gradPhi, source);
    }

    return tPhiSource;
}

// ************************************************************************* //
