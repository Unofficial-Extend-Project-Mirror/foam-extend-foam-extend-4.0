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

#include "rotationalConstraint.H"
#include "sixDOFODE.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rotationalConstraint, 0);
    defineRunTimeSelectionTable(rotationalConstraint, word);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rotationalConstraint::rotationalConstraint
(
    const word& name,
    const dictionary& dict,
    const sixDOFODE& sixDOF
)
:
    name_(name),
    sixDOF_(sixDOF)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rotationalConstraint::~rotationalConstraint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::rotationalConstraint::name() const
{
    return name_;
}


const Foam::sixDOFODE& Foam::rotationalConstraint::sixDOF() const
{
    return sixDOF_;
}


Foam::autoPtr<Foam::rotationalConstraint> Foam::rotationalConstraint::New
(
    const word& name,
    const dictionary& dict,
    const sixDOFODE& sixDOF
)
{
    const word constraintType(dict.lookup("type"));

    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(constraintType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "rotationalConstraint::New"
            "\n("
            "\n    const word& name,"
            "\n    const dictionary& dict,"
            "\n)"
        )   << "Unknown rotation constraint type: " << constraintType
            << endl << endl
            << "Valid rotation constraint types are: " << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<rotationalConstraint>
    (
        cstrIter()
        (
            name,
            dict,
            sixDOF
        )
    );
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const rotationalConstraint& rc
)
{
    os << rc.name_ << nl << token::BEGIN_BLOCK << nl;

    rc.write(os);

    os << token::END_BLOCK << endl;

    os.check("Ostream& operator<<(Ostream&, const rotationalConstraint&");

    return os;
}


// ************************************************************************* //
