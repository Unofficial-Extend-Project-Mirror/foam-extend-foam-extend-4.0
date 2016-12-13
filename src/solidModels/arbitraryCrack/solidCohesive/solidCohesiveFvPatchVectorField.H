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
    solidCohesiveFvPatchVectorField

Description
    Cohesive law fv patch field for arbitrary crack procedure,
    where the mode mixity is allowed to vary depending on how
    the crack opens.

    This procedure allows sigmaMax, tauMax, GIc, and GIIc
    to be independently specified.

    The cohesive law surface field is looked up from rheology.

    The dissipated fracture energy is constantly monitored, by
    integrating in time.

Author
    Philip Cardiff UCD
    based on original Tukovic cohesive law boundary condition.

SourceFiles
    solidCohesiveFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidCohesiveFvPatchVectorField_H
#define solidCohesiveFvPatchVectorField_H

#include "fvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "nonLinearGeometry.H"
#include "cohesiveFvPatch.H"
#include "Switch.H"
#include "tractionBoundaryGradient.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class solidCohesiveFvPatch Declaration
\*---------------------------------------------------------------------------*/

class solidCohesiveFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{

    // Private Data

        //- Name of the displacement field
        const word fieldName_;

        //- Cohesive law - now part of constitutiveModel
        //cohesiveZoneModel* cohesiveZoneModelPtr_;

        //- Under-relaxation factor
        scalar relaxationFactor_;

        //- Current traction
        vectorField traction_;

        //- Minimal unloading relative separation distance
        //scalar minUnloadingSeparationDistance_;

        //- If yes, contact procedure will be used to calculate
        //  cohesive traction for the negative normal separation distance
        Switch contact_;

        //- if a face has cracked
        //Field<bool> cracked_;
        scalarField cracked_;

        //- tractions
        scalarField curTractionN_;
        scalarField oldTractionN_;
        scalarField curTractionS_;
        scalarField oldTractionS_;

        //- deltas
        scalarField deltaN_;
        scalarField oldDeltaN_;
        scalarField deltaS_;
        scalarField oldDeltaS_;
        scalarField unloadingDeltaEff_;

        //- energies
        scalarField currentGI_;
        scalarField oldGI_;
        scalarField currentGII_;
        scalarField oldGII_;


        //- Current time index
        label curTimeIndex_;

        // Penalty factor and returning penalty factor
        // for penalty contact
        scalar* penaltyFactorPtr_;
        scalar penaltyScale_;
        scalar frictionCoeff_; // Coulomb friction coefficient

        // If yes, cohesive traction will be calculated
        // using old separation distance
        Switch explicitSeparationDistance_;
        scalarField curDeltaN_;
        scalarField curDeltaS_;

        // Global patch material field needed for multimaterial cohesive laws
        bool updateGlobalPatchMaterials_;
        scalarField globalPatchMaterials_;

        //- Is it a non linear solver
        nonLinearGeometry::nonLinearType nonLinear_;

        //- if it is an orthropic solver
        Switch orthotropic_;


     // Private member functions

        //- penalty factor
        void calcPenaltyFactor();
        virtual scalar penaltyFactor()
        {
            if (!penaltyFactorPtr_)
            {
                calcPenaltyFactor();
            }

            return *penaltyFactorPtr_;
        }


public:

    //- Runtime type information
    TypeName("solidCohesive");


    // Constructors

        //- Construct from patch and internal field
        solidCohesiveFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidCohesiveFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  solidCohesiveFvPatchVectorField onto a new patch
        solidCohesiveFvPatchVectorField
        (
            const solidCohesiveFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidCohesiveFvPatchVectorField
        (
            const solidCohesiveFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new solidCohesiveFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solidCohesiveFvPatchVectorField
        (
            const solidCohesiveFvPatchVectorField&,
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
                new solidCohesiveFvPatchVectorField(*this, iF)
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
            //const cohesiveZoneModel& law() const;

            //- Return current traction
            // const vectorField& traction() const
            vectorField& traction()
            {
                return traction_;
            }

            //- Normal traction
            //tmp<scalarField> normalTraction() const;

            //- Relative separation distance
            //tmp<scalarField> relativeSeparationDistance() const;

            //- Effective traction
            // tmp<scalarField> effectiveTraction
            // (
            //     const vectorField& traction,
            //     const vectorField& normal
            // ) const;

            //- Effective traction
            // scalar effectiveTraction
            // (
            //     const vector& traction,
            //     const vector& normal
            // ) const;

            //- New faces
            // const scalar& beta() const
            // {
            //     return beta_;
            // }

            //- crack and damage field for visualisation
            //- cracked faces set to 2, damaged but not cracked set to 1
            tmp<scalarField> crackingAndDamage() const;

            //- G field for visualisation
            tmp<scalarField> GI() const;
            tmp<scalarField> GII() const;

            //- print number of damaged and cracked faces
            bool cracking();

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
