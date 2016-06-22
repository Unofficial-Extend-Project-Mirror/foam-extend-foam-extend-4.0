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
    solidCohesiveFixedModeMixFvPatchVectorField

Description
    Cohesive law fv patch field for arbitrary crack procedure,
    where the mode mixity is fixed of a damaged face is fixed
    at the initiation mode mixity.
    Fixing the mode-mixity is only a valid approach when the
    it is clear that the mode-mixity does not change or if the
    damage zones are small.

    This boundary was originally called cohesiveLawFvPatchVectorField.
    There are a couple of changes:
    cohesive law is looked up from constitutive model;
    only Dugdale cohesive law is allowed.

    This method is based on:
    Camacho GT and Ortiz M, Computational modelling of impact damage in
    brittle materials, Int J Solids Structures, 33(20-22):2899-2938, 1996.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    solidCohesiveFixedModeMixFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidCohesiveFixedModeMixFvPatchVectorField_H
#define solidCohesiveFixedModeMixFvPatchVectorField_H

#include "fvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "cohesiveFvPatch.H"
//#include "cohesiveLaw.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class cohesiveLawFvPatch Declaration
\*---------------------------------------------------------------------------*/

class solidCohesiveFixedModeMixFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{

    // Private Data

        //- Name of the displacement field
        const word fieldName_;

        //- Name of rheology model
        //const word rheologyName_;

        //- Cohesive law
        //cohesiveLaw* cohesiveLawPtr_;

        //- Under-relaxation factor
        scalar relaxationFactor_;

        //- Current traction
        vectorField traction_;

        //- Current traction
        vectorField initiationTraction_;

        //- Current separation distance
        vectorField separationDistance_;

        //- Current separation distance
        vectorField oldSeparationDistance_;

        //- Current unloading separation distance
        vectorField unloadingSeparationDistance_;

        //- Minimal unloading relative separation distance
        scalar minUnloadingSeparationDistance_;

        //- Ratio between shear and normal critical traction
        scalar beta_;

        //- If yes, cohesive traction will be calculated
        //  using old separation distance
        Switch explicitSeparationDistance_;

        //- If yes, contact procedure will be used to calculate
        //  cohesive traction for the negative normal separation distance
        Switch contact_;

        //- Constant which determins slope of negative (contact)
        //  part of the traction separation curve
        scalar contactConstant_;

        //- Current time index
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("solidCohesiveFixedModeMix");


    // Constructors

        //- Construct from patch and internal field
        solidCohesiveFixedModeMixFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidCohesiveFixedModeMixFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  solidCohesiveFixedModeMixFvPatchVectorField onto a new patch
        solidCohesiveFixedModeMixFvPatchVectorField
        (
            const solidCohesiveFixedModeMixFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidCohesiveFixedModeMixFvPatchVectorField
        (
            const solidCohesiveFixedModeMixFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new solidCohesiveFixedModeMixFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solidCohesiveFixedModeMixFvPatchVectorField
        (
            const solidCohesiveFixedModeMixFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new solidCohesiveFixedModeMixFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return cohesive law
            //const cohesiveLaw& law() const;

            //- Return current traction
            vectorField& initiationTraction()
            {
                return initiationTraction_;
            }

            //- Return const reference to current traction
            const vectorField& traction() const
            {
                return traction_;
            }

            //- Return reference to current traction
            vectorField& traction()
            {
                return traction_;
            }

            //- Normal traction
            //tmp<scalarField> normalTraction() const;

            //- Relative separation distance
            tmp<scalarField> relativeSeparationDistance() const;

            //- Effective traction
            //tmp<scalarField> effectiveTraction
            //(
            //    const vectorField& traction,
            //    const vectorField& normal
            //) const;

            //- Effective traction
            // scalar effectiveTraction
            // (
            //     const vector& traction,
            //     const vector& normal
            // ) const;

            //- New faces
            const scalar& beta() const
            {
                return beta_;
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
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchVectorField&) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
