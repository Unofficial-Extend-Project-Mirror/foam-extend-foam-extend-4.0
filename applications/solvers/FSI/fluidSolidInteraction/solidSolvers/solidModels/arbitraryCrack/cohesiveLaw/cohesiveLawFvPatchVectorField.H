/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    cohesiveLawFvPatchVectorField

Description
    Cohesive law fv patch field

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    cohesiveLawFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef cohesiveLawFvPatchVectorField_H
#define cohesiveLawFvPatchVectorField_H

#include "fvPatchFields.H"

// #include "directionMixedFvPatchFields.H"
#include "directionMixedDisplacementFvPatchVectorField.H"

#include "cohesiveFvPatch.H"
#include "simpleCohesiveLaw.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class cohesiveLawFvPatch Declaration
\*---------------------------------------------------------------------------*/

class cohesiveLawFvPatchVectorField
:
    public directionMixedDisplacementFvPatchVectorField
{

    // Private Data

        //- Cohesive law
        simpleCohesiveLaw* cohesiveLawPtr_;

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
    TypeName("cohesiveLaw");


    // Constructors

        //- Construct from patch and internal field
        cohesiveLawFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cohesiveLawFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  cohesiveLawFvPatchVectorField onto a new patch
        cohesiveLawFvPatchVectorField
        (
            const cohesiveLawFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        cohesiveLawFvPatchVectorField
        (
            const cohesiveLawFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new cohesiveLawFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        cohesiveLawFvPatchVectorField
        (
            const cohesiveLawFvPatchVectorField&,
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
                new cohesiveLawFvPatchVectorField(*this, iF)
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
            const simpleCohesiveLaw& law() const;

            //- Return current traction
            vectorField& initiationTraction()
            {
                return initiationTraction_;
            }

            //- Return current traction
            const vectorField& traction() const
            {
                return traction_;
            }

            //- Normal traction
            tmp<scalarField> normalTraction() const;

            //- Relative separation distance
            tmp<scalarField> relativeSeparationDistance() const;

            //- Effective traction
            tmp<scalarField> effectiveTraction
            (
                const vectorField& traction,
                const vectorField& normal
            ) const;

            //- Effective traction
            scalar effectiveTraction
            (
                const vector& traction,
                const vector& normal
            ) const;

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
