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
    Foam::overlapGgiEnthalpyJumpFvPatchField

Description
    Foam::overlapGgiEnthalpyJumpFvPatchField

Author
    Ilaria De Dominicis, General Electric Power, (March 2016)

Contributor
    Hrvoje Jasak, Wikki Ltd.
    Gregor Cvijetic, FMENA Zagreb.

GE CONFIDENTIAL INFORMATION 2016 General Electric Company. All Rights Reserved

SourceFiles
    overlapGgiEnthalpyJumpFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef overlapGgiEnthalpyJumpFvPatchField_H
#define overlapGgiEnthalpyJumpFvPatchField_H

#include "jumpOverlapGgiFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class overlapGgiEnthalpyJumpFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class overlapGgiEnthalpyJumpFvPatchField
:
    public jumpOverlapGgiFvPatchField<Type>
{
    // Private data

        //- Is the patch on the rotating side?
        Switch rotating_;

        //- "jump" field
        Field<Type> jump_;


public:

    //- Runtime type information
    TypeName("overlapGgiEnthalpyJump");


// Constructors

        //- Construct from patch and internal field
        overlapGgiEnthalpyJumpFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        overlapGgiEnthalpyJumpFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given overlapGgiEnthalpyJumpFvPatchField
        //  onto a new patch
        overlapGgiEnthalpyJumpFvPatchField
        (
            const overlapGgiEnthalpyJumpFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        overlapGgiEnthalpyJumpFvPatchField
        (
            const overlapGgiEnthalpyJumpFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new overlapGgiEnthalpyJumpFvPatchField<Type>(*this)
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
                new overlapGgiEnthalpyJumpFvPatchField<Type>(*this, iF)
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
void overlapGgiEnthalpyJumpFvPatchField<scalar>::updateCoeffs();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "overlapGgiEnthalpyJumpFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
