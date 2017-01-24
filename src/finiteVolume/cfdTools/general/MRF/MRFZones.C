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

#include "MRFZones.H"
#include "foamTime.H"
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

    forAll (*this, i)
    {
        operator[](i).relativeFlux(MRFZonesPhiCorr);
    }

    return tMRFZonesPhiCorr;
}


Foam::tmp<Foam::surfaceScalarField> Foam::MRFZones::meshPhi() const
{
    tmp<surfaceScalarField> tMRFZonesPhiCorr
    (
        new surfaceScalarField
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
            dimensionedScalar("zero", dimVolume/dimTime, 0)
        )
    );
    surfaceScalarField& MRFZonesPhiCorr = tMRFZonesPhiCorr();

    forAll (*this, i)
    {
        operator[](i).meshPhi(MRFZonesPhiCorr);
    }

    return tMRFZonesPhiCorr;
}


void Foam::MRFZones::addCoriolis(fvVectorMatrix& UEqn) const
{
    forAll (*this, i)
    {
        operator[](i).addCoriolis(UEqn);
    }
}


void Foam::MRFZones::relativeFlux(surfaceScalarField& phi) const
{
    forAll (*this, i)
    {
        operator[](i).relativeFlux(phi);
    }
}


void Foam::MRFZones::absoluteFlux(surfaceScalarField& phi) const
{
    forAll (*this, i)
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
    forAll (*this, i)
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
    forAll (*this, i)
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
    forAll (*this, i)
    {
        operator[](i).absoluteFlux(rho, phi);
    }
}


void Foam::MRFZones::relativeVelocity(volVectorField& U) const
{
    forAll (*this, i)
    {
        operator[](i).relativeVelocity(U);
    }
}


void Foam::MRFZones::absoluteVelocity(volVectorField& U) const
{
    forAll (*this, i)
    {
        operator[](i).absoluteVelocity(U);
    }
}


void Foam::MRFZones::correctBoundaryVelocity(volVectorField& U) const
{
    forAll (*this, i)
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

    // Due to gradient cacheing, must take a tmp field
    // HJ, 22/Apr/2016
    tmp<volVectorField> tgradPhi = fvc::grad(phi);
    const volVectorField& gradPhi = tgradPhi();

    forAll (*this, i)
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

    // Due to gradient cacheing, must take a tmp field
    // HJ, 22/Apr/2016
    tmp<volTensorField> tgradPhi = fvc::grad(phi);
    const volTensorField& gradPhi = tgradPhi();

    forAll (*this, i)
    {
        operator[](i).Su(phi, gradPhi, source);
    }

    return tPhiSource;
}


// ************************************************************************* //
