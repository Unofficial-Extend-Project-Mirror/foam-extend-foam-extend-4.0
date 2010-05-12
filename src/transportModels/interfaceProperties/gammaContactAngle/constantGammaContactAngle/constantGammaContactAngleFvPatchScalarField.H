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
    Foam::constantGammaContactAngleFvPatchScalarField

Description
    A constant gammaContactAngle scalar boundary condition 
    (gammaContactAngleFvPatchScalarField)

SourceFiles
    constantGammaContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef constantGammaContactAngleFvPatchScalarField_H
#define constantGammaContactAngleFvPatchScalarField_H

#include "gammaContactAngleFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class constantGammaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class constantGammaContactAngleFvPatchScalarField
:
    public gammaContactAngleFvPatchScalarField
{
    // Private data

        //- Equilibrium contact angle
        scalar theta0_;


public:

    //- Runtime type information
    TypeName("constantGammaContactAngle");


    // Constructors

        //- Construct from patch and internal field
        constantGammaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        constantGammaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  constantGammaContactAngleFvPatchScalarField
        //  onto a new patch
        constantGammaContactAngleFvPatchScalarField
        (
            const constantGammaContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        constantGammaContactAngleFvPatchScalarField
        (
            const constantGammaContactAngleFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new constantGammaContactAngleFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        constantGammaContactAngleFvPatchScalarField
        (
            const constantGammaContactAngleFvPatchScalarField&,
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
                new constantGammaContactAngleFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Return the equilibrium contact-angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nHat
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
