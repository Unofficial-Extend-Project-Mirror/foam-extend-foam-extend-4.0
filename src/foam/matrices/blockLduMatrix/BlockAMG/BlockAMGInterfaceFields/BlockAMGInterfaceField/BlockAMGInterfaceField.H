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
    Foam::BlockAMGInterfaceField

Description
     Abstract base class for AMG agglomerated interface fields.

Author
    Klas Jareteg, 2013-02-08

SourceFiles
    BlockAMGInterfaceField.C
    newBlockAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef BlockAMGInterfaceField_H
#define BlockAMGInterfaceField_H

#include "BlockLduInterfaceField.H"
#include "AMGInterface.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BlockAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockAMGInterfaceField
:
    public BlockLduInterfaceField<Type>
{
    // Private data

        //- Local reference cast into the interface
        const AMGInterface& interface_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockAMGInterfaceField(const BlockAMGInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const BlockAMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("BlockAMGInterfaceField");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            BlockAMGInterfaceField,
            lduInterface,
            (
                const AMGInterface& AMGCp,
                const BlockLduInterfaceField<Type>& fineInterface
            ),
            (AMGCp, fineInterface)
        );


    // Selectors

        //- Return a pointer to a new interface created on freestore given
        //  the fine interface
        static autoPtr<BlockAMGInterfaceField<Type> > New
        (
            const AMGInterface& AMGCp,
            const BlockLduInterfaceField<Type>& fineInterface
        );


    // Constructors

        //- Construct from AMG interface and fine level interface field
        BlockAMGInterfaceField
        (
            const AMGInterface& AMGCp,
            const BlockLduInterfaceField<Type>&
        )
        :
            BlockLduInterfaceField<Type>(AMGCp),
            interface_(AMGCp)
        {}


    //- Destructor
    virtual ~BlockAMGInterfaceField<Type>()
    {}


    // Member Functions

        // Agglomeration

            // Klas Jareteg: 2013-02-06
            //- Agglomerating for the CoeffField fine-level coefficients
            virtual tmp<CoeffField<Type> > agglomerateBlockCoeffs
            (
                const CoeffField<Type>& fineCoeffs
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockAMGInterfaceField.C"
#   include "newBlockAMGInterfaceField.C"
#endif

#define makeBlockAMGInterfaceField(BlockAMGInterfaceFieldType, typeBlockAMGInterfaceFieldType) \
                                                                              \
defineNamedTemplateTypeNameAndDebug(typeBlockAMGInterfaceFieldType, 0);       \
                                                                              \
addToRunTimeSelectionTable(BlockAMGInterfaceFieldType, typeBlockAMGInterfaceFieldType, lduInterface);

#define makeBlockAMGInterfaceFields(blockAMGInterfaceFieldType)               \
                                                                              \
makeBlockAMGInterfaceField(blockScalarAMGInterfaceField, blockAMGInterfaceFieldType##Scalar); \
makeBlockAMGInterfaceField(blockVectorAMGInterfaceField, blockAMGInterfaceFieldType##Vector); \

//HJ, hacked! Needs terminal specialisation.  HJ, 16/Mar/2016
//  makeBlockAMGInterfaceField(blockTensorAMGInterfaceField, blockAMGInterfaceFieldType##Tensor);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
