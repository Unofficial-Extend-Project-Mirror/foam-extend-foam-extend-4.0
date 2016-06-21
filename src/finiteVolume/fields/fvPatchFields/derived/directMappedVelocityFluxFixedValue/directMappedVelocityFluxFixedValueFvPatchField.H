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
    Foam::directMappedVelocityFluxFixedValueFvPatchField

Description
    Recycles the velocity and flux at a patch to this patch

SourceFiles
    directMappedVelocityFluxFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef directMappedVelocityFluxFixedValueFvPatchField_H
#define directMappedVelocityFluxFixedValueFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "directMappedFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class directMappedVelocityFluxFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

class directMappedVelocityFluxFixedValueFvPatchField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Name of flux field
        word phiName_;

public:

    //- Runtime type information
    TypeName("directMappedVelocityFlux");


    // Constructors

        //- Construct from patch and internal field
        directMappedVelocityFluxFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        directMappedVelocityFluxFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  directMappedVelocityFluxFixedValueFvPatchField
        //  onto a new patch
        directMappedVelocityFluxFixedValueFvPatchField
        (
            const directMappedVelocityFluxFixedValueFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        directMappedVelocityFluxFixedValueFvPatchField
        (
            const directMappedVelocityFluxFixedValueFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new directMappedVelocityFluxFixedValueFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        directMappedVelocityFluxFixedValueFvPatchField
        (
            const directMappedVelocityFluxFixedValueFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new directMappedVelocityFluxFixedValueFvPatchField(*this, iF)
            );
        }


    // Member functions

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
