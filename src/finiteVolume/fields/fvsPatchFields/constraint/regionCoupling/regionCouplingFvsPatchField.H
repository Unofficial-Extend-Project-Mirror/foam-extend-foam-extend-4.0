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
    regionCouplingFvsPatchField

Description
    Region couple patch field

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Note
    The name of the class needs to be different from the name of the patch
    (regionCouple) to avoid constraint patch behaviour

SourceFiles
    regionCouplingFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCouplingFvsPatchField_H
#define regionCouplingFvsPatchField_H

#include "coupledFvsPatchField.H"
#include "regionCoupleLduInterfaceField.H"
#include "regionCoupleFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class regionCouplingFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class regionCouplingFvsPatchField
:
    public coupledFvsPatchField<Type>
{
    // Private data

        //- Local reference cast into the regionCoupling patch
        const regionCoupleFvPatch& regionCouplePatch_;

        //- Name of remote field to couple to
        word remoteFieldName_;

        //- Matrix update buffer
        mutable scalarField matrixUpdateBuffer_;

        //- Original patch field.  Required for correct evaluation
        //  in harmonic averaging
        mutable Field<Type> originalPatchField_;

        //- Current time index used to store originalPatchField_
        mutable label curTimeIndex_;


protected:

        //- Set remote field name
        void setRemoteFieldName(const word& name)
        {
            remoteFieldName_ = name;
        }

        //- Return non-const access to matrixUpdateBuffer
        void setMatrixUpdateBuffer(const scalarField& field) const
        {
            matrixUpdateBuffer_ = field;
        }

        //- Return original patch field
        const Field<Type>& originalPatchField() const;


public:

    //- Runtime type information
    TypeName("regionCoupling");


    // Constructors

        //- Construct from patch and internal field
        regionCouplingFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCouplingFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given regionCouplingFvsPatchField
        //  onto a new patch
        regionCouplingFvsPatchField
        (
            const regionCouplingFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        regionCouplingFvsPatchField
        (
            const regionCouplingFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvsPatchField<Type> > clone() const
        {
            return tmp<fvsPatchField<Type> >
            (
                new regionCouplingFvsPatchField<Type>(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchField<Type> > clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchField<Type> >
            (
                new regionCouplingFvsPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access functions

            //- Return remote field name
            const word& remoteFieldName() const
            {
                return remoteFieldName_;
            }

            //- Return contents of a matrix update buffer
            const scalarField& matrixUpdateBuffer() const
            {
                return matrixUpdateBuffer_;
            }

            //- Return a named shadow patch field
            template<class LookupField, class LookupType>
            const typename LookupField::PatchFieldType&
            lookupShadowPatchField
            (
                const word& name,
                const LookupField* = NULL,
                const LookupType* = NULL
            ) const;

            //- Return a shadow patch
            const regionCoupleFvPatch& regionCouplePatch() const
            {
                return regionCouplePatch_;
            }

            //- Return shadow patch field
            const regionCouplingFvsPatchField<Type>& shadowPatchField() const;


        // Evaluation functions


        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "regionCouplingFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
