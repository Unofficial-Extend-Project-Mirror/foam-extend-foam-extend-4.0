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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Bilinear cohesive law.

\*---------------------------------------------------------------------------*/

#include "bilinearCohesiveLaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bilinearCohesiveLaw, 0);
    addToRunTimeSelectionTable(cohesiveLaw, bilinearCohesiveLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::bilinearCohesiveLaw::bilinearCohesiveLaw
(
    const word& cohesiveLawName,
    const dictionary& dict
)
:
    cohesiveLaw(cohesiveLawName, dict),
    sigma1_(cohesiveLawCoeffs().lookup("sigma1")),
    delta1_(cohesiveLawCoeffs().lookup("delta1")),
    deltaC_((2.0*GIc() - sigmaMax()*delta1())/sigma1())
{}


Foam::bilinearCohesiveLaw::bilinearCohesiveLaw
(
    const bilinearCohesiveLaw& blcl
)
:
    cohesiveLaw(blcl),
    sigma1_(cohesiveLawCoeffs().lookup("sigma1")),
    delta1_(cohesiveLawCoeffs().lookup("delta1")),
    deltaC_(blcl.deltaC_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bilinearCohesiveLaw::~bilinearCohesiveLaw()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return current holding traction
Foam::scalar Foam::bilinearCohesiveLaw::traction(scalar delta) const
{
    if (delta > deltaC().value())
    {
        return 0.0;
    }

    else if (delta < 0)
    {
        return sigmaMax().value();
    }

    else if (delta > delta1().value())
    {
        return sigma1().value()*(1.0 - (delta-delta1().value())/(deltaC().value()-delta1().value()));      
    }

    return sigmaMax().value() + (sigma1().value()-sigmaMax().value())*delta/delta1().value();
}

// ************************************************************************* //
