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

#include "SchaefferFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SchaefferFrictionalStress, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        SchaefferFrictionalStress,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SchaefferFrictionalStress::SchaefferFrictionalStress
(
    const dictionary& dict
)
:
    frictionalStressModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SchaefferFrictionalStress::~SchaefferFrictionalStress()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::SchaefferFrictionalStress::
frictionalPressure
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Fr,
    const dimensionedScalar& eta,
    const dimensionedScalar& p
) const
{
    return
        dimensionedScalar("1e24", dimensionSet(1, -1, -2, 0, 0), 1e24)
       *pow(Foam::max(alpha - alphaMinFriction, scalar(0)), 10.0);
}


Foam::tmp<Foam::volScalarField> Foam::SchaefferFrictionalStress::
frictionalPressurePrime
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Fr,
    const dimensionedScalar& eta,
    const dimensionedScalar& p
) const
{
    return
        dimensionedScalar("1e25", dimensionSet(1, -1, -2, 0, 0), 1e25)
       *pow(Foam::max(alpha - alphaMinFriction, scalar(0)), 9.0);
}


Foam::tmp<Foam::volScalarField> Foam::SchaefferFrictionalStress::muf
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax,
    const volScalarField& pf,
    const volTensorField& D,
    const dimensionedScalar& phi
) const
{
    const scalar I2Dsmall = 1.0e-15;

    // Creating muf assuming it should be 0 on the boundary which may not be
    // true
    tmp<volScalarField> tmuf
    (
        new volScalarField
        (
            IOobject
            (
                "muf",
                alpha.mesh().time().timeName(),
                alpha.mesh()
            ),
            alpha.mesh(),
            dimensionedScalar("muf", dimensionSet(1, -1, -1, 0, 0), 0.0)
        )
    );

    volScalarField& muff = tmuf();

    forAll (D, celli)
    {
        if (alpha[celli] > alphaMax.value()-5e-2)
        {
            muff[celli] =
                0.5*pf[celli]*sin(phi.value())
               /(
                    sqrt(1.0/6.0*(sqr(D[celli].xx() - D[celli].yy())
                  + sqr(D[celli].yy() - D[celli].zz())
                  + sqr(D[celli].zz() - D[celli].xx()))
                  + sqr(D[celli].xy()) + sqr(D[celli].xz())
                  + sqr(D[celli].yz())) + I2Dsmall
                );
        }
    }

    return tmuf;
}


// ************************************************************************* //
