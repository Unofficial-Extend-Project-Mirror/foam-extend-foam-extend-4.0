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

Class
    sixDOFbodies

Description
    6-DOF solver for multiple bodies

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "sixDOFbodies.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sixDOFbodies::setBodies()
{
    // Find if duplicate name existes
    forAll (names_, bodyI)
    {
        for
        (
            label otherBody = bodyI + 1;
            otherBody < names_.size();
            otherBody++
        )
        {
            if (names_[bodyI] == names_[otherBody])
            {
                FatalErrorIn("sixDOFbodies::setBodies()")
                    << "Duplicate names of bodies: this is not allowed"
                    << exit(FatalError);
            }
        }
    }

    odes_.setSize(names_.size());
    solvers_.setSize(names_.size());

    forAll (names_, bodyI)
    {
        odes_.set
        (
            bodyI,
            new sixDOFqODE
            (
                IOobject
                (
                    names_[bodyI],
                    runTime_.timeName(),
                    runTime_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            )
        );

        solvers_.set
        (
            bodyI,
            ODESolver::New
            (
                lookup("solver"),
                odes_[bodyI]
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDOFbodies::sixDOFbodies
(
    const Time& runTime
)
:
    IOdictionary
    (
        IOobject
        (
            "sixDOFsolverDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    runTime_(runTime),
    odes_(),
    solvers_(),
    names_(lookup("bodies"))
{
    setBodies();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sixDOFbodies::solve()
{
    const scalar tol = readScalar(lookup("eps"));

    forAll (odes_, bodyI)
    {
        Info << "Solving 6-DOF for " << names_[bodyI] << " in time"
         << tab << "T = " << runTime_.value() << " s" << endl;

        solvers_[bodyI].solve
        (
            runTime_.value(),
            runTime_.value() + runTime_.deltaT().value(),
            tol,
            runTime_.deltaT().value()
        );
    }
}


const Foam::wordList& Foam::sixDOFbodies::names() const
{
    return names_;
}


const Foam::PtrList<Foam::sixDOFqODE>& Foam::sixDOFbodies::operator()() const
{
    return odes_;
}


// ************************************************************************* //
