/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    Foam::processorBlockGAMGInterfaceField

Description
    GAMG agglomerated processor interface field.

SourceFiles
    processorBlockGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef processorBlockGAMGInterfaceField_H
#define processorBlockGAMGInterfaceField_H

#include "BlockGAMGInterfaceField.H"
#include "processorGAMGInterface.H"
#include "processorBlockLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class processorBlockGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class processorBlockGAMGInterfaceField
:
    public BlockGAMGInterfaceField<Type>,
    public processorBlockLduInterfaceField<Type>
{
    // Private data

        //- Local reference cast into the processor interface
        const processorGAMGInterface& procInterface_;

        //- Is the transform required
        bool doTransform_;

        //- Rank of component for transformation
        int rank_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        processorBlockGAMGInterfaceField
        (
            const processorBlockGAMGInterfaceField&
        );

        //- Disallow default bitwise assignment
        void operator=(const processorBlockGAMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("processor");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        processorBlockGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const BlockLduInterfaceField<Type>& fineInterfaceField
        );


    // Destructor

        virtual ~processorBlockGAMGInterfaceField();


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
                processorBlockLduInterfaceField<Type>::transformCoupleField
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
                processorBlockLduInterfaceField<Type>::transformCoupleField(f);
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

            //- Return rank of component for transform
            virtual int rank() const
            {
                return rank_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "processorBlockGAMGInterfaceField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
