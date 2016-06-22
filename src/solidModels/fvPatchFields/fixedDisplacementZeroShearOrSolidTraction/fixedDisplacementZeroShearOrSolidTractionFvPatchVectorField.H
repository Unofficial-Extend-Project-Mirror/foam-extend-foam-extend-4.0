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
    fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField

Description
    This boundary condition can act either as fixedDisplacementZeroShear
    or solidTraction dictated by the dispTimeSeries.
    For times when dispTimeSeries is >0 then the boundary is
    fixedDisplacementZeroShear.
    For times when dispTimeSeries is 0 then the boundary is
    solidTraction.
    This boundary essentially allows you to change a boundary from
    displacement to traction or vice versa.

SourceFiles
    fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField_H
#define fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField_H

#include "fvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "nonLinearGeometry.H"
#include "Switch.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{
    // Private data

        //- Name of the displacement field
        const word fieldName_;

        //- Traction
        vectorField traction_;

        //- Pressure
        scalarField pressure_;

        //- Displacement
        vectorField displacement_;

        //- Normal to fix for displacement phase
        vector fixedNormal_;

        //- Is it a non linear solver
        nonLinearGeometry::nonLinearType nonLinear_;

        //- Is it an orthropic solver
        Switch orthotropic_;

        //- The time series - 0 for traction, anything else for displacement.
        interpolationTable<scalar> timeSeries_;


    // Private Member Functions


public:

    //- Runtime type information
    TypeName("fixedDisplacementZeroShearOrSolidTraction");


    // Constructors

        //- Construct from patch and internal field
        fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField
        //  onto a new patch
        fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField
        (
            const fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField
                (
                    *this
                    )
            );
        }

        //- Construct as copy setting internal field reference
        fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField
        (
            const fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField&,
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
                new fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField
                (
                    *this,
                    iF
                    )
            );
        }

    // Destructor

        virtual ~fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField()
        {}


    // Member functions

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

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
