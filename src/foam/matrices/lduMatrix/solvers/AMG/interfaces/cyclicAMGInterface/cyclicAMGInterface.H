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
    Foam::cyclicAMGInterface

Description
    AMG agglomerated cyclic interface.

SourceFiles
    cyclicAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMGInterface_H
#define cyclicAMGInterface_H

#include "AMGInterface.H"
#include "cyclicLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cyclicAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class cyclicAMGInterface
:
    public AMGInterface,
    virtual public cyclicLduInterface
{
    // Private data

        //- Reference tor the cyclicLduInterface from which this is
        //  agglomerated
        const cyclicLduInterface& fineCyclicInterface_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cyclicAMGInterface(const cyclicAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const cyclicAMGInterface&);


public:

    //- Runtime type information
    TypeName("cyclic");


    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        cyclicAMGInterface
        (
            const lduPrimitiveMesh& lduMesh,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& localRestrictAddressing,
            const labelField& neighbourRestrictAddressing
        );


    //- Destructor
    virtual ~cyclicAMGInterface();


    // Member Functions

        // Access

            //- Return true if interface is coupled
            virtual bool coupled() const
            {
                return true;
            }


        // Interface transfer functions

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const;

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& iF
            ) const;


        //- Cyclic interface functions

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return fineCyclicInterface_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return fineCyclicInterface_.reverseT();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
