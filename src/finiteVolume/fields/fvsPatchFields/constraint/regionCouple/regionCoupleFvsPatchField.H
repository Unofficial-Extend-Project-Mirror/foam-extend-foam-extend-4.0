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
    regionCoupleFvsPatchField

Description
    Region couple patch field

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    regionCoupleFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupleFvsPatchField_H
#define regionCoupleFvsPatchField_H

#include "coupledFvsPatchField.H"
#include "regionCoupleLduInterfaceField.H"
#include "regionCoupleFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class regionCoupleFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class regionCoupleFvsPatchField
:
    public coupledFvsPatchField<Type>
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
        regionCoupleFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupleFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given regionCoupleFvsPatchField
        //  onto a new patch
        regionCoupleFvsPatchField
        (
            const regionCoupleFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        regionCoupleFvsPatchField
        (
            const regionCoupleFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvsPatchField<Type> > clone() const
        {
            return tmp<fvsPatchField<Type> >
            (
                new regionCoupleFvsPatchField<Type>(*this)
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
                new regionCoupleFvsPatchField<Type>(*this, iF)
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
            const regionCoupleFvsPatchField<Type>& shadowPatchField() const;


        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "regionCoupleFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
