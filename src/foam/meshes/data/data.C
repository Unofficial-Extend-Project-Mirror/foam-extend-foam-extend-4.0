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

#include "data.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::data::debug(Foam::debug::debugSwitchFromDict("data", false));

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::data::data(const objectRegistry& obr)
:
    IOdictionary
    (
        IOobject
        (
            "data",
            obr.time().system(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    prevTimeIndex_(0)
{
    set("solverPerformance", dictionary());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::data::solverPerformanceDict() const
{
    return subDict("solverPerformance");
}


void Foam::data::setSolverPerformance
(
    const word& name,
    const lduSolverPerformance& sp
) const
{
    dictionary& dict = const_cast<dictionary&>(solverPerformanceDict());

    List<lduSolverPerformance> perfs;

    if (prevTimeIndex_ != this->time().timeIndex())
    {
        // Reset solver performance between iterations
        prevTimeIndex_ = this->time().timeIndex();
        dict.clear();
    }
    else
    {
        dict.readIfPresent(name, perfs);
    }

    // Append to list
    perfs.setSize(perfs.size()+1, sp);

    dict.set(name, perfs);
}


void Foam::data::setSolverPerformance
(
    const lduSolverPerformance& sp
) const
{
    setSolverPerformance(sp.fieldName(), sp);
}


// ************************************************************************* //
