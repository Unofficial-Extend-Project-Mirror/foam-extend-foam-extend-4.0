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
    Foam::lduPrimitiveMesh

Description
    Simplest contrete lduMesh which stores the addressing needed by the
    lduMatrix.  Addressing is contained locally

\*---------------------------------------------------------------------------*/

#ifndef lduPrimitiveMesh_H
#define lduPrimitiveMesh_H

#include "lduMesh.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class lduPrimitiveMesh Declaration
\*---------------------------------------------------------------------------*/

class lduPrimitiveMesh
:
    public lduMesh,
    public lduAddressing
{
    // Private data

        //- Lower addressing
        labelList lowerAddr_;

        //- Upper addressing
        labelList upperAddr_;

        //- Patch to internal addressing
        labelListList patchAddr_;

        //- List of pointers for each patch
        //  with only those pointing to interfaces being set
        lduInterfacePtrsList interfaces_;

        //- Patch field evaluation schedule.
        //  Note this does must be held as a copy.  HJ, 20/Feb/2009
        lduSchedule patchSchedule_;

        //- Communicator to use for parallel communication
        const label comm_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        lduPrimitiveMesh(const lduPrimitiveMesh&);

        //- Disallow default bitwise assignment
        void operator=(const lduPrimitiveMesh&);


public:

    // Constructors

        //- Construct from components as copies
        //  Interfaces are added upon construction.  HJ, 20/Feb/2009
        lduPrimitiveMesh
        (
            const label nCells,
            const unallocLabelList& l,
            const unallocLabelList& u,
            const label comm
        )
        :
            lduAddressing(nCells),
            lowerAddr_(l),
            upperAddr_(u),
            comm_(comm)
        {}


        //- Construct from components and re-use storage as specified.
        lduPrimitiveMesh
        (
            const label nCells,
            labelList& l,
            labelList& u,
            const label comm,
            bool reUse
        )
        :
            lduAddressing(nCells),
            lowerAddr_(l, reUse),
            upperAddr_(u, reUse),
            comm_(comm)
        {}


    // Destructor

        virtual ~lduPrimitiveMesh()
        {}


    // Member Functions

        //- Add interfaces.  Constructor helper
        void addInterfaces
        (
            const lduInterfacePtrsList& interfaces,
            const labelListList& pa,
            const lduSchedule& ps
        )
        {
            interfaces_ = interfaces;
            patchAddr_ = pa;
            patchSchedule_ = ps;
        }

        //- Clear interfaces
        void clearInterfaces()
        {
            // Clear interface pointers: they are not deleted
            // by the lduInterfacePtrsList (derived from UPtrList)
            // HJ, 5/Feb/2016
            forAll (interfaces_, i)
            {
                if (interfaces_.set(i))
                {
                    delete interfaces_(i);
                }
            }
        }

        //- Return number of interfaces
        virtual label nPatches() const
        {
            return patchAddr_.size();
        }


        // Access

            //- Return ldu addressing
            virtual const lduAddressing& lduAddr() const
            {
                return *this;
            }

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const
            {
                return interfaces_;
            }

            //- Return Lower addressing
            virtual const unallocLabelList& lowerAddr() const
            {
                return lowerAddr_;
            }

            //- Return Upper addressing
            virtual const unallocLabelList& upperAddr() const
            {
                return upperAddr_;
            }

            //- Return patch addressing
            virtual const unallocLabelList& patchAddr(const label i) const
            {
                return patchAddr_[i];
            }

            //- Return patch evaluation schedule
            virtual const lduSchedule& patchSchedule() const
            {
                return patchSchedule_;
            }

            //- Return communicator used for parallel communication
            virtual label comm() const
            {
                return comm_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
