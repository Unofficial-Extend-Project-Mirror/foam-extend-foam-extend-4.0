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
    cohesiveZoneIncrementalFvPatchVectorField

Description

SourceFiles
    cohesiveZoneIncrementalFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef cohesiveZoneIncrementalFvPatchVectorField_H
#define cohesiveZoneIncrementalFvPatchVectorField_H

#include "fvPatchFields.H"
#include "directionMixedDisplacementFvPatchVectorField.H"
#include "simpleCohesiveLaw.H"
#include "Switch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class cohesiveZoneIncrementalFvPatchVectorField
:
    public directionMixedDisplacementFvPatchVectorField
{
    // Private data

        //- Accumulated displacement
        vectorField totRefValue_;

        //- Name of the total displacement field
        const word totalFieldName_;

        //- Cohesive law
        simpleCohesiveLaw* cohesiveLawPtr_;

        //- Crack indicator
        scalarField crackIndicator_;

        //- CRAZE indicator
        scalarField crazeIndicator_;

        //- Under-relaxation factor
        scalar relaxationFactor_;

        //- Separation distance
        scalarField separationDistance_;

        //- Olde separation distance
        scalarField oldSeparationDistance_;

        //- Current unloading separation distance
        scalarField unloadingSeparationDistance_;

        //- Use explicit spearation distance
        Switch explicitSeparationDistance_;

        //- Current time index
        label curTimeIndex_;

        //- Initiation traction
        vectorField initiationTraction_;

        //- Break only one face per updateCrack
        Switch breakOnlyOneFace_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("cohesiveZoneIncremental");


    // Constructors

        //- Construct from patch and internal field
        cohesiveZoneIncrementalFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cohesiveZoneIncrementalFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  cohesiveZoneIncrementalFvPatchVectorField onto
        //  a new patch
        cohesiveZoneIncrementalFvPatchVectorField
        (
            const cohesiveZoneIncrementalFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new cohesiveZoneIncrementalFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        cohesiveZoneIncrementalFvPatchVectorField
        (
            const cohesiveZoneIncrementalFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new cohesiveZoneIncrementalFvPatchVectorField(*this, iF)
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
            const simpleCohesiveLaw& law() const
            {
                return *cohesiveLawPtr_;
            }

            //- Return total displacement
            const vectorField& totRefValue() const
            {
                return totRefValue_;
            }

            //- Return crack indicator
            const scalarField& crackIndicator() const
            {
                return crackIndicator_;
            }

            //- Return CRAZE indicator
            const scalarField& crazeIndicator() const
            {
                return crazeIndicator_;
            }

            //- Return separation distance
            const scalarField& separationDistance() const
            {
                return separationDistance_;
            }

            //- Return relative separation distance
            const tmp<scalarField> relSeparationDistance() const
            {
                return tmp<scalarField>
                (
                    new scalarField(separationDistance_/law().deltaC().value())
                );
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
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Update crack
            label updateCrack();

//             //- Evaluate the patch field
//             virtual void evaluate
//             (
//                 const Pstream::commsTypes commsType=Pstream::blocking
//             );

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
