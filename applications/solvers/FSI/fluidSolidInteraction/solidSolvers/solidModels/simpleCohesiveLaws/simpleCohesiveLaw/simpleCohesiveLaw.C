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
    Virtual base class for cohesive law.

\*---------------------------------------------------------------------------*/

#include "simpleCohesiveLaw.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleCohesiveLaw, 0);
    defineRunTimeSelectionTable(simpleCohesiveLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::simpleCohesiveLaw> Foam::simpleCohesiveLaw::New
(
    const word& simpleCohesiveLawName,
    const dictionary& dict
)
{
    Info << "Selecting cohesive law: " << simpleCohesiveLawName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(simpleCohesiveLawName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "simpleCohesiveLaw::New(const word& simpleCohesiveLawName, "
            "const dictionary& dict)"
        )   << "Unknown cohesive law " << simpleCohesiveLawName
            << endl << endl
            << "Valid cohesive laws are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<simpleCohesiveLaw>(cstrIter()(simpleCohesiveLawName, dict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::simpleCohesiveLaw::simpleCohesiveLaw
(
    const word& simpleCohesiveLawName,
    const dictionary& dict
)
:
    simpleCohesiveLawCoeffs_(dict.subDict(simpleCohesiveLawName + "Coeffs")),
    GIc_(simpleCohesiveLawCoeffs_.lookup("GIc")),
    sigmaMax_(simpleCohesiveLawCoeffs_.lookup("sigmaMax"))
{}


Foam::simpleCohesiveLaw::simpleCohesiveLaw
(
    const simpleCohesiveLaw& cl
)
:
    refCount(),
    simpleCohesiveLawCoeffs_(cl.simpleCohesiveLawCoeffs_),
    GIc_(cl.GIc_),
    sigmaMax_(cl.sigmaMax_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleCohesiveLaw::~simpleCohesiveLaw()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::simpleCohesiveLaw::writeDict(Ostream& os) const
{
    os.writeKeyword(word(type() + "Coeffs"))
        << simpleCohesiveLawCoeffs();
}


// ************************************************************************* //
