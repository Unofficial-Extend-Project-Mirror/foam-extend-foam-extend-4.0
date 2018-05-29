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
    Foam::nuRatioOmegaFvPatchScalarField

Description
    Calculate turbulent dissipation rate (epsilon) from the provided eddy
    viscosity ratio. Note: derived from inletOutlet so it can readily be used
    for both inlets and outlets.

    Example of the boundary condition specification:
    \verbatim
        inlet
        {
            type        nuRatioOmega;
            nuRatio     10;             // Eddy viscosity ratio
            kName       k;              // Name of the turbulent energy field
                                        // (k by default)
            phiName     phi;            // Name of the flux field
                                        // (phi by default)
            value       uniform 1;      // placeholder
        }
    \endverbatim

SourceFiles
    nuRatioOmegaFvPatchScalarField.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef nuRatioOmegaFvPatchScalarField_H
#define nuRatioOmegaFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class eddyViscosityRatioTurbulenceDissipationRateFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nuRatioOmegaFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private data

        //- Target viscosity ratio
        scalar nuRatio_;

        //- Name of the turbulent energy field
        word kName_;

        //- Name of the flux field
        word phiName_;


public:

    //- Runtime type information
    TypeName("nuRatioOmega");


    // Constructors

        //- Construct from patch and internal field
        nuRatioOmegaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nuRatioOmegaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nuRatioOmegaFvPatchScalarField
        //  onto a new patch
        nuRatioOmegaFvPatchScalarField
        (
            const nuRatioOmegaFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nuRatioOmegaFvPatchScalarField
        (
            const nuRatioOmegaFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nuRatioOmegaFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        nuRatioOmegaFvPatchScalarField
        (
            const nuRatioOmegaFvPatchScalarField&,
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
                new nuRatioOmegaFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
