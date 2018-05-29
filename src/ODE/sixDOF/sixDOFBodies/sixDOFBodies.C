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

Class
    sixDOFBodies

Description
    6-DOF solver for multiple bodies

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "sixDOFBodies.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sixDOFBodies::setBodies()
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
                FatalErrorIn("sixDOFBodies::setBodies()")
                    << "Found duplicate name: " << names_[bodyI]
                    << " for bodies. This is not allowed."
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
            sixDOFODE::New
            (
                IOobject
                (
                    names_[bodyI],
                    runTime_.timeName(),
                    runTime_,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::AUTO_WRITE
                )
            )
        );

        solvers_.set
        (
            bodyI,
            ODESolver::New(lookup("solver"), odes_[bodyI])
        );

        Info<< "Finished creating " << odes_[bodyI].type()
            << " object for body " << names_[bodyI] << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDOFBodies::sixDOFBodies
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
            IOobject::MUST_READ_IF_MODIFIED,
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

void Foam::sixDOFBodies::solve()
{
    const scalar tol = readScalar(lookup("eps"));

    forAll (odes_, bodyI)
    {
        Info << "Solving 6-DOF for " << names_[bodyI] << " in time"
         << tab << "T = " << runTime_.value() << " s" << endl;

        // Note: set external force and moment needed to initialize the state
        // of the sixDOFODE to correctly take into account multiple calls per
        // time step. Using constant force and moment throughout simulation.
        odes_[bodyI].setExternalForceAndMoment
        (
            dimensionedVector(odes_[bodyI].force()),
            dimensionedVector(odes_[bodyI].moment())
        );

        solvers_[bodyI].solve
        (
            runTime_.value() - runTime_.deltaT().value(),
            runTime_.value(),
            tol,
            runTime_.deltaT().value()
        );
    }
}


const Foam::wordList& Foam::sixDOFBodies::names() const
{
    return names_;
}


const Foam::PtrList<Foam::sixDOFODE>& Foam::sixDOFBodies::operator()() const
{
    return odes_;
}


// ************************************************************************* //
