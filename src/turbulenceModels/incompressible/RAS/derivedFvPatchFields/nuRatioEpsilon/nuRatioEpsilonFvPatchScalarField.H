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
    Foam::nuRatioEpsilonFvPatchScalarField

Description
    Calculate turbulent dissipation rate (epsilon) from the provided eddy
    viscosity ratio. Note: derived from inletOutlet so it can readily be used
    for both inlets and outlets.

    Example of the boundary condition specification:
    \verbatim
        inlet
        {
            type        nuRatioEpsilon;
            nuRatio     10;             // Eddy viscosity ratio
            Cmu         0.09;           // Default value 0.09
            kName       k;              // Name of the turbulent energy field
                                        // (k by default)
            phiName     phi;            // Name of the flux field
                                        // (phi by default)
            value       uniform 1;      // placeholder
        }
    \endverbatim

SourceFiles
    nuRatioEpsilonFvPatchScalarField.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef nuRatioEpsilonFvPatchScalarField_H
#define nuRatioEpsilonFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class eddyViscosityRatioTurbulenceDissipationRateFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nuRatioEpsilonFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private data

        //- Target viscosity ratio
        scalar nuRatio_;

        //- Cmu constant
        scalar Cmu_;

        //- Name of the turbulent energy field
        word kName_;

        //- Name of the flux field
        word phiName_;


public:

    //- Runtime type information
    TypeName("nuRatioEpsilon");


    // Constructors

        //- Construct from patch and internal field
        nuRatioEpsilonFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nuRatioEpsilonFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nuRatioEpsilonFvPatchScalarField
        //  onto a new patch
        nuRatioEpsilonFvPatchScalarField
        (
            const nuRatioEpsilonFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nuRatioEpsilonFvPatchScalarField
        (
            const nuRatioEpsilonFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nuRatioEpsilonFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        nuRatioEpsilonFvPatchScalarField
        (
            const nuRatioEpsilonFvPatchScalarField&,
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
                new nuRatioEpsilonFvPatchScalarField
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
