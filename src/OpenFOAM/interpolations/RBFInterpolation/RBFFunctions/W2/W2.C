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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "W2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(W2, 0);
    addToRunTimeSelectionTable(RBFFunction, W2, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::W2::W2(const scalar radius)
:
    RBFFunction(),
    radius_(radius)
{}


Foam::W2::W2(const dictionary& dict)
:
    RBFFunction(),
    radius_(readScalar(dict.lookup("radius")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::W2::~W2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::W2::weights
(
    const vectorField& points,
    const vector& controlPoint
) const
{
    scalarField dist = mag(points - controlPoint);

    scalarField RBF(dist.size());

    RBF = neg(dist - radius_)*
        Foam::max(pow4(1 - (dist/radius_)), scalar(0))*(1 + 4*(dist/radius_));

    return RBF;
}

// ************************************************************************* //
