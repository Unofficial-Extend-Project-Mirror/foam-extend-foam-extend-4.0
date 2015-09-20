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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline const Time& multiSolver::multiDictRegistry() const
{
    return multiDictRegistry_;
}


inline const IOdictionary& multiSolver::multiControlDict() const
{
    return multiControlDict_;
}


inline const wordList& multiSolver::prefixes() const
{
    return prefixes_;
}


inline const label& multiSolver::superLoop() const
{
    return superLoop_;
}


inline const word& multiSolver::currentSolverDomain() const
{
    return currentSolverDomain_;
}


inline const multiSolver::initialStartFromControls&
    multiSolver::initialStartFrom() const
{
    return initialStartFrom_;
}


inline const word& multiSolver::startDomain() const
{
    return startDomain_;
}


inline const scalar& multiSolver::initialStartTime() const
{
    return initialStartTime_;
}


inline const multiSolver::finalStopAtControls& multiSolver::finalStopAt() const
{
    return finalStopAt_;
}


inline const word& multiSolver::endDomain() const
{
    return endDomain_;
}


inline const scalar& multiSolver::finalEndTime() const
{
    return finalEndTime_;
}


inline const label& multiSolver::startSuperLoop() const
{
    return startSuperLoop_;
}


inline const label& multiSolver::endSuperLoop() const
{
    return endSuperLoop_;
}


inline const bool& multiSolver::multiDictsRunTimeModifiable() const
{
    return multiDictsRunTimeModifiable_;
}


inline const scalar& multiSolver::globalTimeOffset() const
{
    return globalTimeOffset_;
}


inline scalar& multiSolver::globalTimeOffset()
{
    return globalTimeOffset_;
}


inline const label& multiSolver::globalIndex() const
{
    return globalIndex_;
}


inline const multiSolver::startFromControls& multiSolver::startFrom() const
{
    return startFrom_;
}

inline const multiSolver::stopAtControls& multiSolver::stopAt() const
{
    return stopAt_;
}


inline const scalar& multiSolver::startTime() const
{
    return startTime_;
}

inline const scalar& multiSolver::endTime() const
{
    return endTime_;
}


inline scalar& multiSolver::endTime()
{
    return endTime_;
}


inline const label& multiSolver::iterations() const
{
    return iterations_;
}


inline const scalar& multiSolver::elapsedTime() const
{
    return elapsedTime_;
}


inline const wordList& multiSolver::storeFields() const
{
    return storeFields_;
}


inline const label& multiSolver::purgeWriteSuperLoops() const
{
    return purgeWriteSuperLoops_;
}


inline const scalar& multiSolver::deltaT() const
{
    return deltaT_;
}


inline scalar& multiSolver::deltaT()
{
    return deltaT_;
}


void Foam::multiSolver::setFinished()
{
    forcedEnd_ = true;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
