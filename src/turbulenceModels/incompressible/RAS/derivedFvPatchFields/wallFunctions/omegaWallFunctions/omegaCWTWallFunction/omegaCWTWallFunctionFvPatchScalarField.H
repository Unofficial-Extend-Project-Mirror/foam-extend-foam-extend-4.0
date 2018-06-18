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
    Foam::incompressible::RASModels::omegaCWTWallFunctionFvPatchScalarField

Description
    Improved boundary condition for specific dissipation, taking into account
    non-equilibrium effects and are sensitive to flow unsteadiness and presence
    of the pressure gradient.

    Reference (bibtex entry):

    @article{popovacHanjalic2007,
        author = {Popovac, M. and Hanjali\'{c}, K.,}
        title = {{Compound Wall Treatment for RANS Computation of Complex
                  Turbulent Flows and Heat Transfer}},
        journal = {Flow Turbulence Combust},
        year = {2007},
        pages = {78--177}
        DOI = {10.1007/s10494-006-9067-x}
    }

SourceFiles
    omegaCWTWallFunctionFvPatchScalarField.C

Author
    Filip Sutalo, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef omegaCWTWallFunctionFvPatchScalarField_H
#define omegaCWTWallFunctionFvPatchScalarField_H

#include "fixedInternalValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
           Class omegaCWTWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class omegaCWTWallFunctionFvPatchScalarField
:
    public fixedInternalValueFvPatchScalarField
{
    // Private data

         //- Name of pressure field
        word pName_;

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

        //- beta1 coefficient
        scalar beta1_;


    // Private member functions

        //- Check the type of the patch
        void checkType();


public:

    //- Runtime type information
    TypeName("omegaCWTWallFunction");


    // Constructors

        //- Construct from patch and internal field
        omegaCWTWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        omegaCWTWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  omegaCWTWallFunctionFvPatchScalarField
        //  onto a new patch
        omegaCWTWallFunctionFvPatchScalarField
        (
            const omegaCWTWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        omegaCWTWallFunctionFvPatchScalarField
        (
            const omegaCWTWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new omegaCWTWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        omegaCWTWallFunctionFvPatchScalarField
        (
            const omegaCWTWallFunctionFvPatchScalarField&,
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
                new omegaCWTWallFunctionFvPatchScalarField(*this, iF)
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
