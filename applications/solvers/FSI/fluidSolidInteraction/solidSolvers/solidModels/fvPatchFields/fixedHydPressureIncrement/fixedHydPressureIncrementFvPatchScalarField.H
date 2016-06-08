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
    fixedHydPressureIncrementFvPatchScalarField

Description
    Fixed hydrostatic pressure boundary condition for the elastic solver

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

SourceFiles
    fixedHydPressureIncrementFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedHydPressureIncrementFvPatchScalarField_H
#define fixedHydPressureIncrementFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class fixedHydPressureIncrementFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedHydPressureIncrementFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

    // Private Data


public:

    //- Runtime type information
    TypeName("fixedHydPressureIncrement");


    // Constructors

        //- Construct from patch and internal field
        fixedHydPressureIncrementFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedHydPressureIncrementFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedHydPressureIncrementFvPatchScalarField onto a new patch
        fixedHydPressureIncrementFvPatchScalarField
        (
            const fixedHydPressureIncrementFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedHydPressureIncrementFvPatchScalarField
        (
            const fixedHydPressureIncrementFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedHydPressureIncrementFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedHydPressureIncrementFvPatchScalarField
        (
            const fixedHydPressureIncrementFvPatchScalarField&,
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
                new fixedHydPressureIncrementFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        // Evaluation functions

            //- Return patch-normal gradient
            virtual tmp<Field<scalar> > snGrad() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<scalar> > gradientBoundaryCoeffs() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
