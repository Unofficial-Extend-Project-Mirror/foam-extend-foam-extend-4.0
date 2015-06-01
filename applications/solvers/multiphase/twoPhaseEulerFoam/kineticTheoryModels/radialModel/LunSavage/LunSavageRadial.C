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

#include "LunSavageRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LunSavageRadial, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        LunSavageRadial,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LunSavageRadial::LunSavageRadial(const dictionary& dict)
:
    radialModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LunSavageRadial::~LunSavageRadial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::LunSavageRadial::g0
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax
) const
{

    return pow(1.0 - alpha/alphaMax, -2.5*alphaMax);
}


Foam::tmp<Foam::volScalarField> Foam::LunSavageRadial::g0prime
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax
) const
{
    return 2.5*alphaMax*alpha*pow(1.0 - alpha, -1.0 - 2.5*alphaMax);
}


// ************************************************************************* //
