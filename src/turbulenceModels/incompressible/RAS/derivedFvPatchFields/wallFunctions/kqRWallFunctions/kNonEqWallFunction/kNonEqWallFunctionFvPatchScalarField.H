/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::incompressible::RASModels::kNonEqWallFunctionFvPatchScalarField

Description
    Boundary condition for turbulence k, Q, and R when using wall functions.
    Acts as a zero gradient condition and takes into account non equilibrium
    effects.

    Reference (bibtex entry):

    @article{sutaloThesis2017,
        author = {\v{S}utalo, F.}
        title = {{Evaluation of Variants of Enhanced Wall Treatment Wall
                  Function in Turbulent Flow Simulations}},
        school = {Faculty of Mechanical Engineering and Naval Architecture,
                  University of Zagreb},
        year = {2017},
    }

SourceFiles
    kNonEqWallFunctionFvPatchScalarField.C

Author
    Filip Sutalo, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef kNonEqWallFunctionFvPatchScalarField_H
#define kNonEqWallFunctionFvPatchScalarField_H

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
          Class kNonEqWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class kNonEqWallFunctionFvPatchScalarField
:
    public zeroGradientFvPatchScalarField
{
    // Private data

        //- Name of velocity field
        word UName_;

        //- Name of turbulence kinetic energy field
        word kName_;

        //- Name of dissipation rate field
        word epsilonName_;

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


    // Private member functions

        //- Check the type of the patch
        void checkType();


public:

    //- Runtime type information
    TypeName("kNonEqWallFunction");


    // Constructors

        //- Construct from patch and internal field
        kNonEqWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        kNonEqWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  kNonEqWallFunctionFvPatchScalarField
        //  onto a new patch
        kNonEqWallFunctionFvPatchScalarField
        (
            const kNonEqWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        kNonEqWallFunctionFvPatchScalarField
        (
            const kNonEqWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new kNonEqWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        kNonEqWallFunctionFvPatchScalarField
        (
            const kNonEqWallFunctionFvPatchScalarField&,
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
                new kNonEqWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

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
