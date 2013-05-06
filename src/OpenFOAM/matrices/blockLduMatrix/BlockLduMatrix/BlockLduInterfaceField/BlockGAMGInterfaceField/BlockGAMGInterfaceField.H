/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::BlockGAMGInterfaceField

Description
     Abstract base class for GAMG agglomerated interface fields.

Author
    Klas Jareteg, 2013-02-08

SourceFiles
    BlockGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef BlockGAMGInterfaceField_H
#define BlockGAMGInterfaceField_H

#include "BlockLduInterfaceField.H"
#include "GAMGInterface.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BlockGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockGAMGInterfaceField
:
    public BlockLduInterfaceField<Type>
{
    // Private data

        //- Local reference cast into the interface
        const GAMGInterface& interface_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockGAMGInterfaceField(const BlockGAMGInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const BlockGAMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("BlockGAMGInterfaceField");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            BlockGAMGInterfaceField,
            lduInterface,
            (
                const GAMGInterface& GAMGCp,
                const BlockLduInterfaceField<Type>& fineInterface
            ),
            (GAMGCp, fineInterface)
        );


    // Selectors

        //- Return a pointer to a new interface created on freestore given
        //  the fine interface
        static autoPtr<BlockGAMGInterfaceField<Type> > New
        (
            const GAMGInterface& GAMGCp,
            const BlockLduInterfaceField<Type>& fineInterface
        );


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        BlockGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const BlockLduInterfaceField<Type>&
        )
        :
            BlockLduInterfaceField<Type>(GAMGCp),
            interface_(GAMGCp)
        {}
            
    // Member Functions

        // Interface matrix update
 
            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField&,
                scalarField&,
                const lduMatrix&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType
            ) const
            {}

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField&,
                scalarField&,
                const lduMatrix&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType
            ) const
            {}


        // Interface block matrix update
 
            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType
            ) const = 0;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockGAMGInterfaceField.C"
#endif

#define makeBlockGAMGInterfaceField(BlockGAMGInterfaceFieldType, typeBlockGAMGInterfaceFieldType)  \
                                                                                                \
defineNamedTemplateTypeNameAndDebug(typeBlockGAMGInterfaceFieldType, 0);                         \
                                                                                                \
addToRunTimeSelectionTable(BlockGAMGInterfaceFieldType, typeBlockGAMGInterfaceFieldType, lduInterface);

#define makeBlockGAMGInterfaceFields(blockGAMGInterfaceFieldType)                                 \
                                                                                                \
makeBlockGAMGInterfaceField(blockScalarGAMGInterfaceField, blockGAMGInterfaceFieldType##Scalar);   \
makeBlockGAMGInterfaceField(blockVectorGAMGInterfaceField, blockGAMGInterfaceFieldType##Vector);   \
makeBlockGAMGInterfaceField(blockTensorGAMGInterfaceField, blockGAMGInterfaceFieldType##Tensor);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
