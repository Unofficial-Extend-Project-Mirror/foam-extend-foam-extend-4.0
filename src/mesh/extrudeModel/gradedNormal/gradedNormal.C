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

#include "gradedNormal.H"
#include "addToRunTimeSelectionTable.H"
#include "BisectionRoot.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace extrudeModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gradedNormal, 0);

addToRunTimeSelectionTable(extrudeModel, gradedNormal, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gradedNormal::gradedNormal(const dictionary& dict)
:
    extrudeModel(typeName, dict),
    thickness_(readScalar(coeffDict_.lookup("thickness"))),
    delta0_(readScalar(coeffDict_.lookup("initialCellLength"))),
    expansionFactor_(1.0)
{
    // Sanity checks
    if (thickness_ <= SMALL)
    {
        FatalErrorIn("gradedNormal(const dictionary&)")
            << "thickness should be positive: " << thickness_
            << exit(FatalError);
    }

    if (delta0_ <= SMALL)
    {
        FatalErrorIn("gradedNormal(const dictionary&)")
            << "initialCellLength should be positive: " << delta0_
            << exit(FatalError);
    }

    const scalar maxExpFactor =
        coeffDict_.lookupOrDefault<scalar>("maxExpansionFactor", 3.0);

    if (maxExpFactor <= SMALL)
    {
        FatalErrorIn("gradedNormal(const dictionary&)")
            << "maxExpansionFactor should be positive: " << maxExpFactor
            << exit(FatalError);
    }

    const scalar bisectionTol =
        coeffDict_.lookupOrDefault<scalar>("bisectionTol", 1e-5);

    if (bisectionTol <= SMALL)
    {
        FatalErrorIn("gradedNormal(const dictionary&)")
            << "bisectionTolerance should be positive: " << bisectionTol
            << exit(FatalError);
    }

    // Create expansion factor equation represented as a function object
    expansionFactorEqn eqn(*this);

    // Calculate the expansionFactor using the bisection algorithm with the
    // default tolerance of 1e-5
    BisectionRoot<expansionFactorEqn> rootFinder
    (
        eqn,
        bisectionTol
    );

    // Search for the root in [0, 3], where upper bound 3 is default
    expansionFactor_ = rootFinder.root
    (
        0.0,
        maxExpFactor
    );

    // Report the result
    Info<< "Calculated expansion factor: " << expansionFactor_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gradedNormal::~gradedNormal()
{}


// * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * * //

point gradedNormal::operator()
(
    const point& surfacePoint,
    const vector& surfaceNormal,
    const label layer
) const
{
    scalar d = 0.0;

    for (label i = 0; i < layer; ++i)
    {
        d += delta0_*pow(expansionFactor_, i + 1);
    }

    return surfacePoint + d*surfaceNormal;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// ************************************************************************* //

