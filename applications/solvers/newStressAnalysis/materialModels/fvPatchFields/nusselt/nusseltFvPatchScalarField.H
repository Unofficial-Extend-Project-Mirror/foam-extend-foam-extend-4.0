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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    nusseltFvPatchScalarField

Description
    Nusselt boundary condition for temperature

SourceFiles
    nusseltFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nusseltFvPatchScalarFields_H
#define nusseltFvPatchScalarFields_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class nusseltFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nusseltFvPatchScalarField
:
    public mixedFvPatchScalarField
{

    // Private data

        //- Name of the displacement field
        const word DTName_;

        //- Reference temperature at infinity
        scalar Tinf_;

        //- Heat transfer coefficient
        scalarField alpha_;


public:

    //- Runtime type information
    TypeName("nusselt");


    // Constructors

        //- Construct from patch and internal field
        nusseltFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nusseltFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nusseltFvPatchScalarField
        //  onto a new patch
        nusseltFvPatchScalarField
        (
            const nusseltFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nusseltFvPatchScalarField
        (
            const nusseltFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nusseltFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nusseltFvPatchScalarField
        (
            const nusseltFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new nusseltFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the temperature at infinity
            scalar Tinf() const
            {
                return Tinf_;
            }

            //- Return reference to the temperature at infinity
            //  to allow adjustment
            scalar& Tinf()
            {
                return Tinf_;
            }

            //- Return the heat transfer coefficient
            const scalarField& alpha() const
            {
                return alpha_;
            }

            //- Return reference to the heat transfer coefficient
            //  to allow adjustment
            scalarField& alpha()
            {
                return alpha_;
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
                const fvPatchField<scalar>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
