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

#include "CloudDistributeTemplate.H"
#include "CloudTemplate.H"
#include "cloud.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::CloudDistribute<ParticleType>::CloudDistribute
(
    const labelList& cellToProc,
    const labelListList& procCellAddressing,
    const labelListList& procFaceAddressing,
    cloud& cloud
)
:
    cloudDistribute(),
    cloud_(reinterpret_cast<Cloud<ParticleType>& >(cloud)),
    transferList_(Pstream::nProcs())
{
    cloud_.split
    (
        cellToProc,
        procCellAddressing,
        procFaceAddressing,
        transferList_
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
void Foam::CloudDistribute<ParticleType>::send
(
    Ostream& toProc,
    const label procIndex
)
{
    // Send transfer list prepared in constructor
    toProc << transferList_[procIndex] << nl;

    // Clear the list for particles to be received
    transferList_[procIndex].clear();
}


template<class ParticleType>
void Foam::CloudDistribute<ParticleType>::receive
(
    Istream& fromProc,
    const label procIndex
)
{
    // Receive particles and transfer into transfer list
    IDLList<ParticleType> newParticles
    (
        fromProc,
        typename ParticleType::iNew(cloud_)
    );

    transferList_[procIndex].DLListBase::transfer(newParticles);
}


template<class ParticleType>
void Foam::CloudDistribute<ParticleType>::rebuild
(
    const PtrList<labelIOList>& cellProcAddressing,
    const PtrList<labelIOList>& faceProcAddressing
)
{
    cloud_.rebuild
    (
        transferList_,
        cellProcAddressing,
        faceProcAddressing
    );
}


// ************************************************************************* //
