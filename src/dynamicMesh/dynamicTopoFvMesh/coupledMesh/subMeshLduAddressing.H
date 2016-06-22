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
    subMeshLduAddressing

Description
    Customized lduAddressing for subMeshes

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef subMeshLduAddressing_H
#define subMeshLduAddressing_H

#include "lduAddressing.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class subMeshLduAddressing Declaration
\*---------------------------------------------------------------------------*/

class subMeshLduAddressing
:
    public lduAddressing
{
    // Private data

        //- Lower as a subList of allOwner
        labelList::subList lowerAddr_;

        //- Upper as a reference to neighbour
        const labelList& upperAddr_;

        //- Number of patches
        label nPatches_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        subMeshLduAddressing(const subMeshLduAddressing&);

        //- Disallow default bitwise assignment
        void operator=(const subMeshLduAddressing&);

public:

    // Constructors

        //- Construct from components
        subMeshLduAddressing(const fvMesh& mesh)
        :
            lduAddressing(mesh.nCells()),
            lowerAddr_
            (
                labelList::subList
                (
                    mesh.faceOwner(),
                    mesh.nInternalFaces()
                )
            ),
            upperAddr_(mesh.faceNeighbour()),
            nPatches_(mesh.boundary().size())
        {}

    // Destructor

        virtual ~subMeshLduAddressing()
        {}

    // Member Functions

        //- Return number of interfaces
        virtual label nPatches() const
        {
            return nPatches_;
        }

        //- Return lower addressing (i.e. lower label = upper triangle)
        virtual const unallocLabelList& lowerAddr() const
        {
            return lowerAddr_;
        }

        //- Return upper addressing (i.e. upper label)
        virtual const unallocLabelList& upperAddr() const
        {
            return upperAddr_;
        }

        //- Return patch addressing
        virtual const unallocLabelList& patchAddr(const label i) const
        {
            FatalErrorIn
            (
                "void subMeshLduAddressing::patchAddr(const label i)"
            )
                << " Illegal request. "
                << abort(FatalError);

            return unallocLabelList::null();
        }

        // Return patch field evaluation schedule
        virtual const lduSchedule& patchSchedule() const
        {
            FatalErrorIn
            (
                "void subMeshLduAddressing::patchSchedule()"
            )
                << " Illegal request. "
                << abort(FatalError);

            return lduSchedule::null();
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
