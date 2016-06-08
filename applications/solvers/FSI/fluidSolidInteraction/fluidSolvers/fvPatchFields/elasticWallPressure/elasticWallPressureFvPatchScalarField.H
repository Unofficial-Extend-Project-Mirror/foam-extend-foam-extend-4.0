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
    elasticWallPressureFvPatchScalarField

Description
    elasticWallPressureFvPatchScalarField

SourceFiles
    elasticWallPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef elasticWallPressureFvPatchScalarField_H
#define elasticWallPressureFvPatchScalarField_H

#include "robinFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class elasticWallPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class elasticWallPressureFvPatchScalarField
:
    public robinFvPatchField<scalar>
{
    // Private data

        //- Local time index
        label timeIndex_;

        //- Pressure applied at the solid side of the interface
        //  in the previous fsi iteration
        scalarField prevPressure_;

        //- Acceleration of the solid side of the interface
        //  in the previous fsi iteration
        vectorField prevAcceleration_;

        //- Old face centres of patch.
        //  Needed for (second order) time integration
        vectorField Fc_;
        vectorField oldFc_;
        vectorField oldOldFc_;

        //- Face centre velocity.
        //  Needed for (second order) time integration
        vectorField U_;
        vectorField oldU_;
        vectorField oldOldU_;

public:

    //- Runtime type information
    TypeName("elasticWallPressure");


    // Constructors

        //- Construct from patch and internal field
        elasticWallPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        elasticWallPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given elasticWallPressureFvPatchScalarField
        //  onto a new patch
        elasticWallPressureFvPatchScalarField
        (
            const elasticWallPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        elasticWallPressureFvPatchScalarField
        (
            const elasticWallPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new elasticWallPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        elasticWallPressureFvPatchScalarField
        (
            const elasticWallPressureFvPatchScalarField&,
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
                new elasticWallPressureFvPatchScalarField(*this, iF)
            );
        }



    // Member functions

        // Access functions

            //- Get previous pressure patch field
            const scalarField& prevPressure() const
            {
                return prevPressure_;
            }


            //- Get previous pressure patch field
            scalarField& prevPressure()
            {
                return prevPressure_;
            }

            //- Get previous pressure patch field
            const vectorField& prevAcceleration() const
            {
                return prevAcceleration_;
            }


            //- Get previous pressure patch field
            vectorField& prevAcceleration()
            {
                return prevAcceleration_;
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

//             //- Return gradient at boundary
//             virtual tmp<Field<scalar> > snGrad() const
//             {
//                 return tmp<Field<scalar> >
//                 (
//                     new scalarField(this->patch().size(), 0)
//                 );
//             }

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Calculate patch flux
            virtual void patchFlux
            (
                GeometricField<scalar, fvsPatchField, surfaceMesh>& flux,
                const fvMatrix<scalar>& matrix
            ) const;

//             //- Evaluate the patch field
//             virtual void evaluate
//             (
//                 const Pstream::commsTypes commsType=Pstream::blocking
//             );

//             //- Return the matrix diagonal coefficients corresponding to the
//             //  evaluation of the value of this patchField with given weights
//             virtual tmp<Field<scalar> > valueInternalCoeffs
//             (
//                 const tmp<scalarField>&
//             ) const;

//             //- Return the matrix source coefficients corresponding to the
//             //  evaluation of the value of this patchField with given weights
//             virtual tmp<Field<scalar> > valueBoundaryCoeffs
//             (
//                 const tmp<scalarField>&
//             ) const;

//             //- Return the matrix diagonal coefficients corresponding to the
//             //  evaluation of the gradient of this patchField
//             virtual tmp<Field<scalar> > gradientInternalCoeffs() const;

//             //- Return the matrix source coefficients corresponding to the
//             //  evaluation of the gradient of this patchField
//             virtual tmp<Field<scalar> > gradientBoundaryCoeffs() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
