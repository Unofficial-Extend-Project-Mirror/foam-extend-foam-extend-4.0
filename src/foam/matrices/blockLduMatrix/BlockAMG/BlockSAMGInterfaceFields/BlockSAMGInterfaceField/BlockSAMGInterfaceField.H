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
    Foam::BlockSAMGInterfaceField

Description
     Abstract base class for AMG selected interface fields.

Author
    Hrvoje Jasak.  All rights reserved

SourceFiles
    BlockSAMGInterfaceField.C
    newBlockSAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef BlockSAMGInterfaceField_H
#define BlockSAMGInterfaceField_H

#include "BlockLduInterfaceField.H"
#include "SAMGInterface.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BlockSAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockSAMGInterfaceField
:
    public BlockLduInterfaceField<Type>
{
    // Private data

        //- Local reference cast into the interface
        const SAMGInterface& interface_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockSAMGInterfaceField(const BlockSAMGInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const BlockSAMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("BlockSAMGInterfaceField");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            BlockSAMGInterfaceField,
            lduInterface,
            (
                const SAMGInterface& SAMGCp,
                const BlockLduInterfaceField<Type>& fineInterface
            ),
            (SAMGCp, fineInterface)
        );


    // Selectors

        //- Return a pointer to a new interface created on freestore given
        //  the fine interface
        static autoPtr<BlockSAMGInterfaceField<Type> > New
        (
            const SAMGInterface& SAMGCp,
            const BlockLduInterfaceField<Type>& fineInterface
        );


    // Constructors

        //- Construct from AMG interface and fine level interface field
        BlockSAMGInterfaceField
        (
            const SAMGInterface& SAMGCp,
            const BlockLduInterfaceField<Type>&
        )
        :
            BlockLduInterfaceField<Type>(SAMGCp),
            interface_(SAMGCp)
        {}


    //- Destructor
    virtual ~BlockSAMGInterfaceField<Type>()
    {}


    // Member Functions

        //- Select the CoeffField fine-level coefficients
        //  for the coarse level
        virtual tmp<CoeffField<Type> > selectBlockCoeffs
        (
            const CoeffField<Type>& fineCoeffs
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockSAMGInterfaceField.C"
#   include "newBlockSAMGInterfaceField.C"
#endif

#define makeBlockSAMGInterfaceField(BlockSAMGInterfaceFieldType, typeBlockSAMGInterfaceFieldType) \
                                                                               \
defineNamedTemplateTypeNameAndDebug(typeBlockSAMGInterfaceFieldType, 0);       \
                                                                               \
addToRunTimeSelectionTable(BlockSAMGInterfaceFieldType, typeBlockSAMGInterfaceFieldType, lduInterface);

#define makeBlockSAMGInterfaceFields(blockSAMGInterfaceFieldType)              \
                                                                               \
makeBlockSAMGInterfaceField(blockScalarSAMGInterfaceField, blockSAMGInterfaceFieldType##Scalar); \
makeBlockSAMGInterfaceField(blockVectorSAMGInterfaceField, blockSAMGInterfaceFieldType##Vector); \

//HJ, hacked! Needs terminal specialisation.  HJ, 16/Mar/2016
//  makeBlockSAMGInterfaceField(blockTensorSAMGInterfaceField, blockSAMGInterfaceFieldType##Tensor);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
