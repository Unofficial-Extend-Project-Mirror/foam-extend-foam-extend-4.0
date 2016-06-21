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

#include "error.H"

#include "RanzMarshall.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(RanzMarshall, 0);

addToRunTimeSelectionTable
(
    heatTransferModel,
    RanzMarshall,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
RanzMarshall::RanzMarshall
(
    const dictionary& dict
)
:
    heatTransferModel(dict),
    heatDict_(dict.subDict(typeName + "Coeffs")),
    preRePrFactor_(readScalar(heatDict_.lookup("preRePrFactor"))),
    ReExponent_(readScalar(heatDict_.lookup("ReExponent"))),
    PrExponent_(readScalar(heatDict_.lookup("PrExponent")))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

RanzMarshall::~RanzMarshall()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool RanzMarshall::heatTransfer() const
{
    return true;
}

scalar RanzMarshall::Nu
(
    const scalar ReynoldsNumber,
    const scalar PrandtlNumber
) const
{
    return 2.0 + preRePrFactor_ * pow(ReynoldsNumber, ReExponent_) * pow(PrandtlNumber, PrExponent_);
}

scalar RanzMarshall::relaxationTime
(
    const scalar liquidDensity,
    const scalar diameter,
    const scalar liquidcL,
    const scalar kappa,
    const scalar ReynoldsNumber,
    const scalar PrandtlNumber
) const
{
    scalar time = liquidDensity*pow(diameter, 2.0)*liquidcL/(6.0*kappa*Nu(ReynoldsNumber, PrandtlNumber));

    time = max(SMALL, time);

    return time;
}

scalar RanzMarshall::fCorrection(const scalar z) const
{
    scalar correct;
    if (z > 0.01)
    {
        if (z < 1.0e+5)
        {
            correct = z/(exp(z) - 1.0);
        }
        else
        {
            correct = SMALL;
        }

    }
    else
    {
        // taylor-expansion of exp(z)...
        correct = 1.0/(1+0.5*z);
    }

    return correct;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
