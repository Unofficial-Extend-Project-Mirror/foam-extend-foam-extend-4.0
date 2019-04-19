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
    Foam::GGIBlockSAMGInterfaceField

Description
    AMG selected GGI interface field.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    GGIBlockSAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef GGIBlockSAMGInterfaceField_H
#define GGIBlockSAMGInterfaceField_H

#include "BlockSAMGInterfaceField.H"
#include "ggiSAMGInterface.H"
#include "GGIBlockLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class GGIBlockSAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class GGIBlockSAMGInterfaceField
:
    public BlockSAMGInterfaceField<Type>,
    public GGIBlockLduInterfaceField<Type>
{
    // Private data

        //- Local reference cast into the ggi interface
        const ggiSAMGInterface& ggiInterface_;

        //- Is the transform required
        bool doTransform_;

        //- Field transfer buffer
        mutable Field<Type> fieldTransferBuffer_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        GGIBlockSAMGInterfaceField(const GGIBlockSAMGInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const GGIBlockSAMGInterfaceField&);


        //- Select type
        template<class Type2>
        void selectBlockType
        (
            Field<Type2>& coarseCoeffs,
            const Foam::Field<Type2>& fineCoeffs
        ) const;


public:

    //- Runtime type information
    TypeName("ggi");


    // Constructors

        //- Construct from SAMG interface and fine level interface field
        GGIBlockSAMGInterfaceField
        (
            const SAMGInterface& SAMGCp,
            const BlockLduInterfaceField<Type>& fineInterfaceField
        );


    //- Destructor
    virtual ~GGIBlockSAMGInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return ggiInterface_.size();
            }


        // Agglomeration

            // Klas Jareteg: 2013-02-06. Moved by HJ, 16/Mar/2016
            //- Agglomerating for the CoeffField fine-level coefficients
            virtual tmp<CoeffField<Type> > selectBlockCoeffs
            (
                const CoeffField<Type>& fineCoeffs
            ) const;


        // Block coupled interface matrix update

            //- Transform given patch component field
            virtual void transformCoupleField
            (
                scalarField& f,
                const direction cmpt
            ) const
            {
                GGIBlockLduInterfaceField<Type>::transformCoupleField
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
                GGIBlockLduInterfaceField<Type>::transformCoupleField(f);
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


        //- Ggi interface functions

            //- Does the interface field perform the transfromation
            virtual bool doTransform() const
            {
                return doTransform_;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return ggiInterface_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return ggiInterface_.reverseT();
            }


        // Transfer buffer access

            //- Return contents of the field transfer buffer
            const Field<Type>& fieldTransferBuffer() const
            {
                return fieldTransferBuffer_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "GGIBlockSAMGInterfaceField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
