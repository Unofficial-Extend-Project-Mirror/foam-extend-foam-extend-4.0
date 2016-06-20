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

#include "tolerances.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tolerances::tolerances(const Time& t, const fileName& dictName)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            t.system(),
            t,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    relaxationFactors_(ITstream("relaxationFactors", tokenList())()),
    solverTolerances_(ITstream("solverTolerances", tokenList())()),
    solverRelativeTolerances_
    (
        ITstream("solverRelativeTolerances", tokenList())()
    )
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool tolerances::read()
{
    if (regIOobject::read())
    {
        const word toleranceSetName(lookup("toleranceSet"));
        const dictionary& toleranceSet(subDict(toleranceSetName));

        if (toleranceSet.found("relaxationFactors"))
        {
            relaxationFactors_ = toleranceSet.subDict("relaxationFactors");
        }

        if (toleranceSet.found("solverTolerances"))
        {
            solverTolerances_ = toleranceSet.subDict("solverTolerances");
        }

        if (toleranceSet.found("solverRelativeTolerances"))
        {
            solverRelativeTolerances_ =
                toleranceSet.subDict("solverRelativeTolerances");
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool tolerances::relax(const word& name) const
{
    return relaxationFactors_.found(name);
}

scalar tolerances::relaxationFactor(const word& name) const
{
    return readScalar(relaxationFactors_.lookup(name));
}

scalar tolerances::solverTolerance(const word& name) const
{
    return readScalar(solverTolerances_.lookup(name));
}

bool tolerances::solverRelativeTolerances() const
{
    return solverRelativeTolerances_.size();
}

scalar tolerances::solverRelativeTolerance(const word& name) const
{
    return readScalar(solverRelativeTolerances_.lookup(name));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
