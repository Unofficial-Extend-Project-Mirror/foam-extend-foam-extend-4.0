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
    Foam::processorSAMGInterface

Description
    AMG selected processor interface.

SourceFiles
    processorSAMGInterface.C
    processorSAMGInterfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef processorSAMGInterface_H
#define processorSAMGInterface_H

#include "SAMGInterface.H"
#include "processorLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class processorSAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class processorSAMGInterface
:
    public SAMGInterface,
    public processorLduInterface
{
    // Private data

        //- Reference tor the processorLduInterface from which this is
        //  agglomerated
        const processorLduInterface& fineProcInterface_;

        //- Communicator to use for parallel communication
        const label comm_;

        //- Message tag used for sending
        const int tag_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        processorSAMGInterface(const processorSAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const processorSAMGInterface&);


public:

    //- Runtime type information
    TypeName("processor");

    // Constructors

        //- Construct from fine-level interface,
        //  local and neighbour row label
        processorSAMGInterface
        (
            const lduPrimitiveMesh& lduMesh,
            const crMatrix& interfaceProlongation,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const crMatrix& nbrInterfaceProlongation
        );


    //- Destructor
    virtual ~processorSAMGInterface();


    // Member Functions

        // Access

            //- Return interface size
            virtual label interfaceSize() const
            {
                return SAMGInterface::size();
            }

            //- Return true if interface is coupled
            virtual bool coupled() const
            {
                return true;
            }

            //- Does this side own the patch ?
            virtual bool master() const
            {
                return myProcNo() < neighbProcNo();
            }


        // Communications support

            //- Return communicator used for parallel communication
            virtual int comm() const
            {
                return comm_;
            }

            //- Return message tag used for sending
            virtual int tag() const
            {
                return tag_;
            }


        // Interface transfer functions

            //- Initialise interface data transfer
            virtual void initTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const;

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const;

            //- Initialise neighbour field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& iF
            ) const;

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& iF
            ) const;


            //- Initialise transfer of prolongation matrix to the interface
            virtual void initProlongationTransfer
            (
                const Pstream::commsTypes commsType,
                const crMatrix& filteredP
            ) const;

            //- Transfer and return prolongation matrix adjacent to
            //  the interface
            virtual autoPtr<crMatrix> prolongationTransfer
            (
                const Pstream::commsTypes commsType,
                const crMatrix& filteredP
            ) const;


        //- Processor interface functions

            //- Return processor number
            virtual int myProcNo() const
            {
                return fineProcInterface_.myProcNo();
            }

            //- Return neigbour processor number
            virtual int neighbProcNo() const
            {
                return fineProcInterface_.neighbProcNo();
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return fineProcInterface_.forwardT();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
