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
    Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a grey-diffuse condition for radiative
    heat flux, \c Qr, for use with the view factor model

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        Qro          | external radiative heat flux | yes    |
        emissivityMode | emissivity mode: solidThermo or lookup | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            greyDiffusiveRadiationViewFactor;
        Qro             uniform 0;
        emissivityMode  solidThermo;
        value           uniform 0;
    }
    \endverbatim

SeeAlso
    Foam::radiation::radiationModel
    Foam::radiation::viewFactor
    Foam::radiationCoupledBase
    Foam::fixedValueFvPatchField

SourceFiles
    greyDiffusiveViewFactorFixedValueFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef greyDiffusiveViewFactorFixedValueFvPatchScalarField_H
#define greyDiffusiveViewFactorFixedValueFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
/*---------------------------------------------------------------------------*\
      Class greyDiffusiveViewFactorFixedValueFvPatchScalarField declaration
\*---------------------------------------------------------------------------*/

class greyDiffusiveViewFactorFixedValueFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- External radiative heat flux
        scalarField Qro_;

        scalarField emissivity_;


public:

    //- Runtime type information
    TypeName("greyDiffusiveRadiationViewFactor");


    // Constructors

        //- Construct from patch and internal field
        greyDiffusiveViewFactorFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        greyDiffusiveViewFactorFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        // greyDiffusiveViewFactorFixedValueFvPatchScalarField onto a new patch
        greyDiffusiveViewFactorFixedValueFvPatchScalarField
        (
            const greyDiffusiveViewFactorFixedValueFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        greyDiffusiveViewFactorFixedValueFvPatchScalarField
        (
            const greyDiffusiveViewFactorFixedValueFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new greyDiffusiveViewFactorFixedValueFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        greyDiffusiveViewFactorFixedValueFvPatchScalarField
        (
            const greyDiffusiveViewFactorFixedValueFvPatchScalarField&,
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
                new greyDiffusiveViewFactorFixedValueFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return external radiative heat flux
            const scalarList& Qro()
            {
                return Qro_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Calculate corresponding emissivity field
            scalarField emissivity() const
            {
                return emissivity_;
            }


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
