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

#include "PstreamGlobals.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Outstanding non-blocking operations.
//! \cond fileScope
DynamicList<MPI_Request> PstreamGlobals::outstandingRequests_;
//! \endcond

// Max outstanding message tag operations.
//! \cond fileScope
int PstreamGlobals::nTags_ = 0;
//! \endcond

// Free'd message tags
//! \cond fileScope
DynamicList<int> PstreamGlobals::freedTags_;
//! \endcond


// Allocated communicators.
//! \cond fileScope
DynamicList<MPI_Comm> PstreamGlobals::MPICommunicators_;
DynamicList<MPI_Group> PstreamGlobals::MPIGroups_;
//! \endcond

void PstreamGlobals::checkCommunicator
(
    const label comm,
    const label otherProcNo
)
{
    if
    (
        comm < 0
     || comm >= PstreamGlobals::MPICommunicators_.size()
    )
    {
        FatalErrorIn
        (
            "PstreamGlobals::checkCommunicator(const label, const label)"
        )   << "otherProcNo:" << otherProcNo << " : illegal communicator "
            << comm << endl
            << "Communicator should be within range 0.."
            << PstreamGlobals::MPICommunicators_.size() - 1
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
