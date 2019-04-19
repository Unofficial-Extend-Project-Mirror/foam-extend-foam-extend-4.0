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
    Foam::incompressible::RASModels::nutMEWTWallFunctionFvPatchScalarField

Description
    Improved boundary condition for turbulent (kinematic) viscosity, taking into
    account convection and pressure gradient terms.

    Reference (bibtex entry):

    @article{sutaloThesis2017,
        author = {\v{S}utalo, F.}
        title = {{Evaluation of Variants of Enhanced Wall Treatment Wall
                  Function in Turbulent Flow Simulations}},
        school = {Faculty of Mechanical Engineering and Naval Architecture,
                  University of Zagreb},
        year = {2017},
    }

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nutCWTWallFunction;
    }
    \endverbatim

See also
    Foam::nutkWallFunctionFvPatchScalarField

SourceFiles
    nutMEWTWallFunctionFvPatchScalarField.C

Author
    Filip Sutalo, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef nutMEWTWallFunctionFvPatchScalarField_H
#define nutMEWTWallFunctionFvPatchScalarField_H

#include "nutkWallFunctionFvPatchScalarField.H"
#include "tolerancesSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
             Class nutMEWTWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutMEWTWallFunctionFvPatchScalarField
:
    public nutkWallFunctionFvPatchScalarField
{
    // Private Data Members

        //- Static tolerance for the A dimensionless parameter
        static const debug::tolerancesSwitch dimlessAFactorTol_;

        //- Name of velocity field
        word UName_;

        //- Name of pressure field
        word pName_;

        //- Name of eddy viscosity field
        word nutName_;


protected:

    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutMEWTWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutMEWTWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutMEWTWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nutMEWTWallFunctionFvPatchScalarField
        //  onto a new patch
        nutMEWTWallFunctionFvPatchScalarField
        (
            const nutMEWTWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutMEWTWallFunctionFvPatchScalarField
        (
            const nutMEWTWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutMEWTWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutMEWTWallFunctionFvPatchScalarField
        (
            const nutMEWTWallFunctionFvPatchScalarField&,
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
                new nutMEWTWallFunctionFvPatchScalarField(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
