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
    Foam::passiveProcessorPolyPatch

Description
    Passive processor polyPatch.

    Passive processor polyPatch is created in load balancing when new processor
    faces are built on each parallel decomposition.  It contains the data
    which is needed to rebuild correct processor boundaries after load balancing

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    passiveProcessorPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef passiveProcessorPolyPatch_H
#define passiveProcessorPolyPatch_H

#include "polyBoundaryMesh.H"
#include "labelField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class passiveProcessorPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class passiveProcessorPolyPatch
:
    public polyPatch
{
    // Private data

        //- My processor number
        int myProcNo_;

        //- Neighbour processor number
        int neighbProcNo_;

        //- Global processor face index
        labelField globalFaceIndex_;


public:

    //- Runtime type information
    TypeName("passiveProcessor");


    // Constructors

        //- Construct from components
        passiveProcessorPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const int myProcNo,
            const int neighbProcNo,
            const labelList& globalFaceIndex
        );

        //- Construct from dictionary
        passiveProcessorPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&
        );

        //- Construct as given the original patch and resetting the
        //  face list and boundary mesh information
        passiveProcessorPolyPatch
        (
            const passiveProcessorPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct as copy
        passiveProcessorPolyPatch(const passiveProcessorPolyPatch&);

        //- Construct as copy, resetting the boundary mesh
        passiveProcessorPolyPatch
        (
            const passiveProcessorPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new passiveProcessorPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new passiveProcessorPolyPatch
                (
                    refCast<const passiveProcessorPolyPatch>(*this),
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~passiveProcessorPolyPatch();


    // Member functions

        //- Return processor number
        int myProcNo() const
        {
            return myProcNo_;
        }

        //- Return neigbour processor number
        int neighbProcNo() const
        {
            return neighbProcNo_;
        }

        //- Return global processor face index
        const labelField& globalFaceIndex() const
        {
            return globalFaceIndex_;
        }

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
