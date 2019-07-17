/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "porousZones.H"
#include "foamTime.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<porousZone>, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousZones::porousZones
(
    const fvMesh& mesh
)
:
    IOPtrList<porousZone>
    (
        IOobject
        (
            "porousZones",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        porousZone::iNew(mesh)
    ),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousZones::addResistance(fvVectorMatrix& UEqn) const
{
    forAll(*this, i)
    {
        operator[](i).addResistance(UEqn);
    }
}


void Foam::porousZones::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    // addResistance for each zone, delaying the correction of the
    // processor BCs of AU
    forAll(*this, i)
    {
        operator[](i).addResistance(UEqn, AU, false);
    }

    // Correct the boundary conditions of the tensorial diagonal to ensure
    // processor bounaries are correctly handled when AU^-1 is interpolated
    // for the pressure equation.
    AU.correctBoundaryConditions();
}


void Foam::porousZones::addHeatResistance
(
    fvScalarMatrix& hTEqn,
    const volScalarField& T,
    volScalarField& Taux,
    volScalarField& Qaux,
    const volVectorField& U,
    const volScalarField& Macro,
    const volScalarField& posFlux
) const
{
    forAll(*this, i)
    {
        operator[](i).addHeatResistance
        (
            hTEqn,
            T,
            Taux,
            Qaux,
            U,
            Macro,
            posFlux
        );
    }
}

// Order cells for Dual Stream model
void Foam::porousZones::macroCellOrder
(
    volScalarField& Taux,
    volScalarField& Macro,
    volScalarField& posFlux,
    const surfaceScalarField& phi
) const
{
    forAll(*this, i)
    {
        operator[](i).macroCellOrder(Taux, Macro, posFlux, phi);
    }
}

bool Foam::porousZones::readData(Istream& is)
{
    clear();

    IOPtrList<porousZone> newLst
    (
        IOobject
        (
            "porousZones",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false     // Don't re-register new zones with objectRegistry
        ),
        porousZone::iNew(mesh_)
    );

    transfer(newLst);

    return is.good();
}


bool Foam::porousZones::writeData(Ostream& os, bool subDict) const
{
    // Write size of list
    os << nl << size();

    // Write beginning of contents
    os << nl << token::BEGIN_LIST;

    // Write list contents
    forAll(*this, i)
    {
        os << nl;
        operator[](i).writeDict(os, subDict);
    }

    // Write end of contents
    os << token::END_LIST << nl;

    // Check state of IOstream
    return os.good();
}


// ************************************************************************* //
