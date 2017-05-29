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
    Foam::incompressible::RASModels::nutCWTWallFunctionFvPatchScalarField

Description
    Improved boundary condition for turbulent (kinematic) viscosity, taking into
    account non-equilibrium effects and are sensitive to flow unsteadiness and
    presence of the pressure gradient.

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
    nutCWTWallFunctionFvPatchScalarField.C

Author
    Filip Sutalo, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef nutCWTWallFunctionFvPatchScalarField_H
#define nutCWTWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
             Class nutCWTWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutCWTWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
protected:

    // Protected data

        //- Name of velocity field
        word UName_;

        //- Name of pressure field
        word pName_;

        //- Name of eddy viscosity field
        word nutName_;

        //- Name of omega field
        word omegaName_;

        //- Cmu coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;

        //- Y+ at the edge of the laminar sublayer
        scalar yPlusLam_;


    // Protected member functions

        //- Check the type of the patch
        virtual void checkType();

        //- Calculate the Y+ at the edge of the laminar sublayer
        virtual scalar calcYPlusLam(const scalar kappa, const scalar E) const;

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutCWTWallTreatment");


    // Constructors

        //- Construct from patch and internal field
        nutCWTWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutCWTWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutCWTWallFunctionFvPatchScalarField
        //  onto a new patch
        nutCWTWallFunctionFvPatchScalarField
        (
            const nutCWTWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutCWTWallFunctionFvPatchScalarField
        (
            const nutCWTWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutCWTWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutCWTWallFunctionFvPatchScalarField
        (
            const nutCWTWallFunctionFvPatchScalarField&,
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
                new nutCWTWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
