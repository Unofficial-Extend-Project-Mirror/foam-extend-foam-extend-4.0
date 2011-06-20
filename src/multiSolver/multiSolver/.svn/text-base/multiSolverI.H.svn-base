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

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline const Time& multiSolver::multiDictRegistry() const
{
    return multiDictRegistry_;
}


inline const IOdictionary multiSolver::multiControlDict() const
{
    return multiControlDict_;
}


inline const label& multiSolver::superLoop() const
{
    return superLoop_;
}


inline const word& multiSolver::currentSolverDomain() const
{
    return currentSolverDomain_;
}


inline const multiSolver::initialStartFromControls& multiSolver::initialStartFrom() const
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


void Foam::multiSolver::setFinished()
{
    finished_ = true;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
