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
    Foam::ProcessorTopology

Description
    Determines processor-processor connection. After instantiation contains
    on all processors the processor-processor connection table.

    *this[procI] gives the list of neighbouring processors.

SourceFiles
    ProcessorTopologyTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef ProcessorTopologyTemplate_H
#define ProcessorTopologyTemplate_H

#include "labelList.H"
#include "lduSchedule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ProcessorTopology Declaration
\*---------------------------------------------------------------------------*/

template<class Patch, class ProcPatch>
class ProcessorTopology
:
    public labelListList
{

private:

    // Private data

        //- Local map from neighbour proc to patchI. Different per processor!
        //  -1 or patchI for connection to procID
        labelList procPatchMap_;

        //- Order in which the patches should be initialised/evaluated
        //  corresponding to the schedule
        lduSchedule patchSchedule_;


    // Private Member Functions

        //- Return all neighbouring processors of this processor. Set
        //  procPatchMap_.
        labelList procNeighbours(const PtrList<Patch>&);

public:

    // Constructors

        //- Construct from boundaryMesh
        ProcessorTopology(const PtrList<Patch>& patches);


    // Member Functions

        //- from neighbour processor to index in boundaryMesh. Local information
        //  (so not same over all processors)
        const labelList& procPatchMap() const
        {
            return procPatchMap_;
        }

        //- Order in which the patches should be initialised/evaluated
        //  corresponding to the schedule
        const lduSchedule& patchSchedule() const
        {
            return patchSchedule_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ProcessorTopologyTemplate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
