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
    Foam::mixingPlaneEnthalpyJumpFvPatchField

Description
    Foam::mixingPlaneEnthalpyJumpFvPatchField

Author
    Ilaria De Dominicis, General Electric Power, (March 2016)

Contributor
    Hrvoje Jasak, Wikki Ltd.

GE CONFIDENTIAL INFORMATION 2016 General Electric Company. All Rights Reserved

SourceFiles
    mixingPlaneEnthalpyJumpFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mixingPlaneEnthalpyJumpFvPatchField_H
#define mixingPlaneEnthalpyJumpFvPatchField_H

#include "jumpMixingPlaneFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class mixingPlaneEnthalpyJumpFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mixingPlaneEnthalpyJumpFvPatchField
:
    public jumpMixingPlaneFvPatchField<Type>
{
    // Private data

        //- Is the patch on the rotating side?
        Switch rotating_;

        //- "jump" field
        Field<Type> jump_;


public:

    //- Runtime type information
    TypeName("mixingPlaneEnthalpyJump");


// Constructors

        //- Construct from patch and internal field
        mixingPlaneEnthalpyJumpFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mixingPlaneEnthalpyJumpFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixingPlaneEnthalpyJumpFvPatchField
        //  onto a new patch
        mixingPlaneEnthalpyJumpFvPatchField
        (
            const mixingPlaneEnthalpyJumpFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        mixingPlaneEnthalpyJumpFvPatchField
        (
            const mixingPlaneEnthalpyJumpFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new mixingPlaneEnthalpyJumpFvPatchField<Type>(*this)
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
                new mixingPlaneEnthalpyJumpFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Is the patch rotating
            bool rotating() const
            {
                return rotating_;
            }

            //- Return the "jump" across the patch as a "half" field
            virtual tmp<Field<Type> > jump() const
            {
                return jump_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


//- Specialisation of the jump-condition for the enthalpy
template<>
void mixingPlaneEnthalpyJumpFvPatchField<scalar>::updateCoeffs();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mixingPlaneEnthalpyJumpFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
