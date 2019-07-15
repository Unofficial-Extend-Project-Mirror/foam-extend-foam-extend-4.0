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
    Foam::incompressible::RASModels::omegaWallFunctionFvPatchScalarField

Description
    Provides a wall function boundary condition/constraint on omega

    Computed value is:

        omega = sqrt(omega_vis^2 + omega_log^2)

    where
        omega_vis = omega in viscous region
        omega_log = omega in logarithmic region

    Model described by Eq.(15) of:
    @verbatim
        Menter, F., Esch, T.
        "Elements of Industrial Heat Transfer Prediction"
        16th Brazilian Congress of Mechanical Engineering (COBEM),
        Nov. 2001
    @endverbatim

Description
    Boundary condition for omega when using wall functions
    on immersed boundary patches

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    immersedBoundaryOmegaWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryOmegaWallFunctionFvPatchScalarField_H
#define immersedBoundaryOmegaWallFunctionFvPatchScalarField_H

#include "omegaWallFunctionFvPatchScalarField.H"
#include "immersedBoundaryFieldBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
     Class immersedBoundaryOmegaWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryOmegaWallFunctionFvPatchScalarField
:
    public omegaWallFunctionFvPatchScalarField,
    public immersedBoundaryFieldBase<scalar>
{
public:

    //- Runtime type information
    TypeName("immersedBoundaryOmegaWallFunction");


    // Constructors

        //- Construct from patch and internal field
        immersedBoundaryOmegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedBoundaryOmegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // immersedBoundaryOmegaWallFunctionFvPatchScalarField
        //  onto a new patch
        immersedBoundaryOmegaWallFunctionFvPatchScalarField
        (
            const immersedBoundaryOmegaWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        immersedBoundaryOmegaWallFunctionFvPatchScalarField
        (
            const immersedBoundaryOmegaWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new immersedBoundaryOmegaWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        immersedBoundaryOmegaWallFunctionFvPatchScalarField
        (
            const immersedBoundaryOmegaWallFunctionFvPatchScalarField&,
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
                new immersedBoundaryOmegaWallFunctionFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


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
                const fvPatchScalarField&,
                const labelList&
            );

            //- Update on mesh motion
            virtual void updateOnMotion();


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patchField
            virtual void evaluate
            (
                const Pstream::commsTypes  = Pstream::blocking
            );

            //- Manipulate a matrix
            virtual void manipulateMatrix(fvScalarMatrix& matrix);


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
