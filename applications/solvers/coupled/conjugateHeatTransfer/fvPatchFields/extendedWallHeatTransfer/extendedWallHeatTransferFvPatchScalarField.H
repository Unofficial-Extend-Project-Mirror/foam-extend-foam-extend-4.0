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
    Foam::extendedWallHeatTransferFvPatchScalarField

Description
    wall heat transfer boundary condition taking into account

    - heat resistance through a solid planar wall

      q = hc*(Tb - Two)

    - heat convection

      q = alpha*(Two - Tinf)

    - grey body radiation to infinite medium

      q = epsilon*sigma*(Two^4 - Tinf^4)

SourceFiles
    extendedWallHeatTransferFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef extendedWallHeatTransferFvPatchScalarField_H
#define extendedWallHeatTransferFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "externalRadiationSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class wallHeatTransferFvPatch Declaration
\*---------------------------------------------------------------------------*/

class extendedWallHeatTransferFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Tinf [K]
        scalar Tinf_;

        //- thermal conductance [W/m^2/K]
        scalar hc_;

        //- heat transfer coefficient [W/m^2/K]
        scalar alpha_;

        //- list of radiative sources
        PtrList<externalRadiationSource> radSources_;

        //- Name of thermal conductivity field
        word KName_;

        //- Take radiation into account for this patch
        bool radiation_;


public:

    //- Runtime type information
    TypeName("extendedWallHeatTransfer");


    // Constructors

        //- Construct from patch and internal field
        extendedWallHeatTransferFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        extendedWallHeatTransferFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given extendedWallHeatTransferFvPatchScalarField
        // onto a new patch
        extendedWallHeatTransferFvPatchScalarField
        (
            const extendedWallHeatTransferFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        extendedWallHeatTransferFvPatchScalarField
        (
            const extendedWallHeatTransferFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new extendedWallHeatTransferFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        extendedWallHeatTransferFvPatchScalarField
        (
            const extendedWallHeatTransferFvPatchScalarField&,
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
                new extendedWallHeatTransferFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return radiation
            bool radiation() const
            {
                return radiation_;
            }

        // Mapping functions

        // Evaluation functions

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
