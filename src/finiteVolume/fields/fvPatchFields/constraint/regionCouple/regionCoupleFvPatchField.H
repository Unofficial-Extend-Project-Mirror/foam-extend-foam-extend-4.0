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
    regionCoupleFvPatchField

Description
    Region couple patch field

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    regionCoupleFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupleFvPatchField_H
#define regionCoupleFvPatchField_H

#include "coupledFvPatchField.H"
#include "regionCoupleLduInterfaceField.H"
#include "regionCoupleFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class regionCoupleFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class regionCoupleFvPatchField
:
    public regionCoupleLduInterfaceField,
    public coupledFvPatchField<Type>
{
    // Private data

        //- Local reference cast into the regionCouple patch
        const regionCoupleFvPatch& regionCouplePatch_;

        //- Name of remote field to couple to
        word remoteFieldName_;

        //- Matrix update buffer
        mutable scalarField matrixUpdateBuffer_;


protected:

        //- Return contents of a matrix update buffer
        const scalarField& matrixUpdateBuffer() const
        {
            return matrixUpdateBuffer_;
        }


public:

    //- Runtime type information
//     TypeName(regionCoupleFvPatch::typeName_());
    TypeName("regionCoupling");


    // Constructors

        //- Construct from patch and internal field
        regionCoupleFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupleFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given regionCoupleFvPatchField onto a new patch
        regionCoupleFvPatchField
        (
            const regionCoupleFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        regionCoupleFvPatchField
        (
            const regionCoupleFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new regionCoupleFvPatchField<Type>(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new regionCoupleFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Return remote field name
            const word& remoteFieldName() const
            {
                return remoteFieldName_;
            }

            //- Return shadow patch field
            const regionCoupleFvPatchField<Type>& shadowPatchField() const;

            //- Return neighbour field given internal cell data
            virtual tmp<Field<Type> > patchNeighbourField() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType
            );

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;


        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "regionCoupleFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
