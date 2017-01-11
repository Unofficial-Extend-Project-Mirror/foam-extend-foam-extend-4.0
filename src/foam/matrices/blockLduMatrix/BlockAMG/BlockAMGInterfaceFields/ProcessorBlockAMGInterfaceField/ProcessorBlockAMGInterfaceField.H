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
    Foam::ProcessorBlockAMGInterfaceField

Description
    AMG agglomerated processor interface field.

SourceFiles
    ProcessorBlockAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef ProcessorBlockAMGInterfaceField_H
#define ProcessorBlockAMGInterfaceField_H

#include "BlockAMGInterfaceField.H"
#include "processorAMGInterface.H"
#include "ProcessorBlockLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class ProcessorBlockAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class ProcessorBlockAMGInterfaceField
:
    public BlockAMGInterfaceField<Type>,
    public ProcessorBlockLduInterfaceField<Type>
{
    // Private data

        //- Local reference cast into the processor interface
        const processorAMGInterface& procInterface_;

        //- Is the transform required
        bool doTransform_;


        // Sending and receiving

            //- Outstanding request
            mutable label outstandingSendRequest_;

            //- Outstanding request
            mutable label outstandingRecvRequest_;

            //- Send buffer.
            mutable Field<Type> sendBuf_;

            //- Receive buffer.
            mutable Field<Type> receiveBuf_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ProcessorBlockAMGInterfaceField
        (
            const ProcessorBlockAMGInterfaceField&
        );

        //- Disallow default bitwise assignment
        void operator=(const ProcessorBlockAMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("processor");


    // Constructors

        //- Construct from AMG interface and fine level interface field
        ProcessorBlockAMGInterfaceField
        (
            const AMGInterface& AMGCp,
            const BlockLduInterfaceField<Type>& fineInterfaceField
        );


    //- Destructor
    virtual ~ProcessorBlockAMGInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return procInterface_.size();
            }


        // Block coupled interface matrix update

            //- Transform given patch component field
            virtual void transformCoupleField
            (
                scalarField& f,
                const direction cmpt
            ) const
            {
                ProcessorBlockLduInterfaceField<Type>::transformCoupleField
                (
                    f,
                    cmpt
                );
            }

            //- Transform neighbour field
            virtual void transformCoupleField
            (
                Field<Type>& f
            ) const
            {
                ProcessorBlockLduInterfaceField<Type>::transformCoupleField(f);
            }

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;


        //- Processor interface functions

            //- Return communicator used for comms
            virtual label comm() const
            {
                return procInterface_.comm();
            }

            //- Return processor number
            virtual int myProcNo() const
            {
                return procInterface_.myProcNo();
            }

            //- Return neigbour processor number
            virtual int neighbProcNo() const
            {
                return procInterface_.neighbProcNo();
            }

            //- Does the interface field perform the transfromation
            virtual bool doTransform() const
            {
                return doTransform_;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return procInterface_.forwardT();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ProcessorBlockAMGInterfaceField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
