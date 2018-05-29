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

#include "solution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class FieldType>
void Foam::solution::cachePrintMessage
(
    const char* message,
    const word& name,
    const FieldType& vf
)
{
    if (solution::debug)
    {
        Info<< "Cache: " << message << token::SPACE << name
            << ", originating from " << vf.name()
            << " event No. " << vf.eventNo()
            << endl;
    }
}


template<class Type>
void Foam::solution::setSolverPerformance
(
    const word& name,
    const BlockSolverPerformance<Type>& sp
) const
{
    List<BlockSolverPerformance<Type> > perfs;

    if (prevTimeIndex_ != this->time().timeIndex())
    {
        // Reset solver performance between iterations
        prevTimeIndex_ = this->time().timeIndex();
        solverPerformance_.clear();
    }
    else
    {
        solverPerformance_.readIfPresent(name, perfs);
    }

    // If storeAllResiduals_ is true, we are storing residual of every iteration
    // inside a single time step. Otherwise, only the first iteration residual
    // and the current iteration residual are required, so the current
    // iteration residual replaces the previous one and only the first iteration
    // residual is always present, VS 2018-02-11
    if (storeAllResiduals_ || perfs.size() < 2)
    {
        // Append to list
        perfs.setSize(perfs.size() + 1, sp);
    }
    else
    {
        perfs.last() = sp;
    }

    solverPerformance_.set(name, perfs);
}


template<class Type>
void Foam::solution::setSolverPerformance
(
    const BlockSolverPerformance<Type>& sp
) const
{
    setSolverPerformance(sp.fieldName(), sp);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
