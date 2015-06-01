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

#include "LunPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LunPressure, 0);

    addToRunTimeSelectionTable
    (
        granularPressureModel,
        LunPressure,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LunPressure::LunPressure(const dictionary& dict)
:
    granularPressureModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LunPressure::~LunPressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::LunPressure::granularPressureCoeff
(
    const volScalarField& alpha,
    const volScalarField& g0,
    const dimensionedScalar& rhoa,
    const dimensionedScalar& e
) const
{

    return rhoa*alpha*(1.0 + 2.0*(1.0 + e)*alpha*g0);
}


Foam::tmp<Foam::volScalarField> Foam::LunPressure::granularPressureCoeffPrime
(
    const volScalarField& alpha,
    const volScalarField& g0,
    const volScalarField& g0prime,
    const dimensionedScalar& rhoa,
    const dimensionedScalar& e
) const
{
    return rhoa*(1.0 + alpha*(1.0 + e)*(4.0*g0 + 2.0*g0prime*alpha));
}

// ************************************************************************* //
