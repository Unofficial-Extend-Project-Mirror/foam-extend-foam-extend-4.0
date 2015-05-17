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

Description
    Virtual base class for cohesive law.

\*---------------------------------------------------------------------------*/

#include "cohesiveLaw.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cohesiveLaw, 0);
    defineRunTimeSelectionTable(cohesiveLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::cohesiveLaw> Foam::cohesiveLaw::New
(
    const word& cohesiveLawName,
    const dictionary& dict
)
{
    Info << "Selecting cohesive law: " << cohesiveLawName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(cohesiveLawName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "cohesiveLaw::New(const word& cohesiveLawName, "
            "const dictionary& dict)"
        )   << "Unknown cohesive law " << cohesiveLawName
            << endl << endl
            << "Valid cohesive laws are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<cohesiveLaw>(cstrIter()(cohesiveLawName, dict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::cohesiveLaw::cohesiveLaw
(
    const word& cohesiveLawName,
    const dictionary& dict
)
:
    cohesiveLawCoeffs_(dict.subDict(cohesiveLawName + "Coeffs")),
    GIc_(cohesiveLawCoeffs_.lookup("GIc")),
    sigmaMax_(cohesiveLawCoeffs_.lookup("sigmaMax"))
{}


Foam::cohesiveLaw::cohesiveLaw
(
    const cohesiveLaw& cl
)
:
    refCount(),
    cohesiveLawCoeffs_(cl.cohesiveLawCoeffs_),
    GIc_(cl.GIc_),
    sigmaMax_(cl.sigmaMax_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cohesiveLaw::~cohesiveLaw()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cohesiveLaw::writeDict(Ostream& os) const
{
    os.writeKeyword(word(type() + "Coeffs"))
        << cohesiveLawCoeffs();
}


// ************************************************************************* //
