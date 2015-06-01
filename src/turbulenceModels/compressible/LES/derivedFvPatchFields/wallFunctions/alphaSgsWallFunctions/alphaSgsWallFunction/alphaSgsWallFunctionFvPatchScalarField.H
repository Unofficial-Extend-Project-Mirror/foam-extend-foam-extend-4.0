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
    alphaSgsWallFunctionFvPatchScalarField

Description
    Boundary condition for thermal diffusivity when using wall functions
    - replicates OpenFOAM v1.5 (and earlier) behaviour

SourceFiles
    alphaSgsWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphaSgsWallFunctionFvPatchScalarField_H
#define alphaSgsWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
           Class alphaSgsWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphaSgsWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Turbulent Prandtl number
        scalar Prt_;


    // Private member functions

        //- Check the type of the patch
        void checkType();


public:

    //- Runtime type information
    TypeName("alphaSgsWallFunction");


    // Constructors

        //- Construct from patch and internal field
        alphaSgsWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphaSgsWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given alphaSgsWallFunctionFvPatchScalarField
        //  onto a new patch
        alphaSgsWallFunctionFvPatchScalarField
        (
            const alphaSgsWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        alphaSgsWallFunctionFvPatchScalarField
        (
            const alphaSgsWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new alphaSgsWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        alphaSgsWallFunctionFvPatchScalarField
        (
            const alphaSgsWallFunctionFvPatchScalarField&,
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
                new alphaSgsWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Evaluate the patchField
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::Pstream::blocking
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
