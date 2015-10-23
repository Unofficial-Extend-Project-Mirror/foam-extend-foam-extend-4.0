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

#include "vanDriestDelta.H"
#include "LESModel.H"
#include "wallDistData.H"
#include "wallPointYPlus.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(vanDriestDelta, 0);
addToRunTimeSelectionTable(LESdelta, vanDriestDelta, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void vanDriestDelta::calcDelta()
{
    const LESModel& lesModel = mesh_.lookupObject<LESModel>("LESProperties");

    const volVectorField& U = lesModel.U();
    const volScalarField& nu = lesModel.nu();
    tmp<volScalarField> nuSgs = lesModel.nuSgs();

    volScalarField ystar
    (
        IOobject
        (
            "ystar",
            mesh_.time().constant(),
            U.db()
        ),
        mesh_,
        dimensionedScalar("ystar", dimLength, GREAT)
    );

    const fvPatchList& patches = mesh_.boundary();
    forAll(patches, patchi)
    {
        if (patches[patchi].isWall())
        {
            const fvPatchVectorField& Uw = U.boundaryField()[patchi];
            const scalarField& nuw = nu.boundaryField()[patchi];
            const scalarField& nuSgsw = nuSgs().boundaryField()[patchi];

            ystar.boundaryField()[patchi] =
                nuw/sqrt((nuw + nuSgsw)*mag(Uw.snGrad()) + VSMALL);
        }
    }

    wallPointYPlus::yPlusCutOff = 500;
    wallDistData<wallPointYPlus> y(mesh_, ystar);

    delta_ = min
    (
        static_cast<const volScalarField&>(geometricDelta_()),
        (kappa_/Cdelta_)*((scalar(1) + SMALL) - exp(-y/ystar/Aplus_))*y
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vanDriestDelta::vanDriestDelta
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dd
)
:
    LESdelta(name, mesh),
    geometricDelta_
    (
        LESdelta::New("geometricDelta", mesh, dd.subDict(type() + "Coeffs"))
    ),
    kappa_(dd.lookupOrDefault<scalar>("kappa", 0.41)),
    Aplus_
    (
        dd.subDict(type() + "Coeffs").lookupOrDefault<scalar>("Aplus", 26.0)
    ),
    Cdelta_
    (
        dd.subDict(type() + "Coeffs").lookupOrDefault<scalar>("Cdelta", 0.158)
    ),
    calcInterval_
    (
        dd.subDict(type() + "Coeffs").lookupOrDefault<label>("calcInterval", 1)
    )
{
    delta_ = geometricDelta_();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void vanDriestDelta::read(const dictionary& d)
{
    const dictionary& dd(d.subDict(type() + "Coeffs"));

    geometricDelta_().read(dd);
    d.readIfPresent<scalar>("kappa", kappa_);
    dd.readIfPresent<scalar>("Aplus", Aplus_);
    dd.readIfPresent<scalar>("Cdelta", Cdelta_);
    dd.readIfPresent<label>("calcInterval", calcInterval_);
    calcDelta();
}


void vanDriestDelta::correct()
{
    if (mesh().time().timeIndex() % calcInterval_ == 0)
    {
        geometricDelta_().correct();
        calcDelta();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
