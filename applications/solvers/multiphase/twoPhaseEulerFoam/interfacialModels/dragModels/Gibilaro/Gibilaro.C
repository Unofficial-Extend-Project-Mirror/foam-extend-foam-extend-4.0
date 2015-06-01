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

#include "Gibilaro.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Gibilaro, 0);

    addToRunTimeSelectionTable
    (
        dragModel,
        Gibilaro,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Gibilaro::Gibilaro
(
    const dictionary& interfaceDict,
    const volScalarField& alpha,
    const phaseModel& phasea,
    const phaseModel& phaseb
)
:
    dragModel(interfaceDict, alpha, phasea, phaseb)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Gibilaro::~Gibilaro()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::Gibilaro::K
(
    const volScalarField& Ur
) const
{
    volScalarField beta = max(scalar(1) - alpha_, scalar(1.0e-6));
    volScalarField bp = pow(beta, -2.8);
    volScalarField Re = max(beta*Ur*phasea_.d()/phaseb_.nu(), scalar(1.0e-3));

    return (17.3/Re + scalar(0.336))*phaseb_.rho()*Ur*bp/phasea_.d();
}


// ************************************************************************* //
