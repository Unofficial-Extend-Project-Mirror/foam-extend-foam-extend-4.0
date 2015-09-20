/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    Foam::incompressible::RASModels::
    immersedBoundaryEpsilonWallFunctionFvPatchScalarField

Description
    Boundary condition for epsilon when using wall functions
    - calculates y+, G, tangential velocity, nut and k
    - each of calculated values is handled separately by the appropriate
      boundary condition
    - epsilon values added directly into the field to act as a constraint

SourceFiles
    immersedBoundaryEpsilonWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryEpsilonWallFunctionFvPatchScalarField_H
#define immersedBoundaryEpsilonWallFunctionFvPatchScalarField_H

#include "immersedBoundaryWallFunctionFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
     Class immersedBoundaryEpsilonWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryEpsilonWallFunctionFvPatchScalarField
:
    public immersedBoundaryWallFunctionFvPatchScalarField
{
    // Private data

        //- Name of velocity field
        word UName_;

        //- Name of turbulence kinetic energy field
        word kName_;

        //- Name of turbulence generation field
        word GName_;

        //- Name of laminar viscosity field
        word nuName_;

        //- Name of turbulent viscosity field
        word nutName_;

        //- Cmu coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;


public:

    //- Runtime type information
    TypeName("immersedBoundaryEpsilonWallFunction");


    // Constructors

        //- Construct from patch and internal field
        immersedBoundaryEpsilonWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedBoundaryEpsilonWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  immersedBoundaryEpsilonWallFunctionFvPatchScalarField
        //  onto a new patch
        immersedBoundaryEpsilonWallFunctionFvPatchScalarField
        (
            const immersedBoundaryEpsilonWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        immersedBoundaryEpsilonWallFunctionFvPatchScalarField
        (
            const immersedBoundaryEpsilonWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new immersedBoundaryEpsilonWallFunctionFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        immersedBoundaryEpsilonWallFunctionFvPatchScalarField
        (
            const immersedBoundaryEpsilonWallFunctionFvPatchScalarField&,
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
                new immersedBoundaryEpsilonWallFunctionFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    //- Destructor
    virtual ~immersedBoundaryEpsilonWallFunctionFvPatchScalarField()
    {}


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patchField
            virtual void evaluate
            (
                const Pstream::commsTypes  = Pstream::blocking
            );


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
