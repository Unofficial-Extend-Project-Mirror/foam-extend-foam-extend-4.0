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

#include "laminarFlameSpeed.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::laminarFlameSpeed> Foam::laminarFlameSpeed::New
(
    const hhuCombustionThermo& ct
)
{
    IOdictionary laminarFlameSpeedDict
    (
        IOobject
        (
            "combustionProperties",
            ct.T().time().constant(),
            ct.T().db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word laminarFlameSpeedType
    (
        laminarFlameSpeedDict.lookup("laminarFlameSpeedCorrelation")
    );

    Info<< "Selecting laminar flame speed correlation "
        << laminarFlameSpeedType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(laminarFlameSpeedType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "laminarFlameSpeed::New(const hhuCombustionThermo&)",
            laminarFlameSpeedDict
        )   << "Unknown laminarFlameSpeed type "
            << laminarFlameSpeedType << endl << endl
            << "Valid laminarFlameSpeed types are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<laminarFlameSpeed>(cstrIter()(laminarFlameSpeedDict, ct));
}


// ************************************************************************* //
