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
    Foam::compressible::RASModels::
        mutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField

Description
    This boundary condition is the same as mutURoughWallFunction. We will simply
    derive from mutURoughWallFunctionFvPatchScalarField to enable backward
    compatibility.

SourceFiles
    mutSpalartAllamarasStandardWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleMutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField_H
#define compressibleMutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField_H

#include "mutURoughWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
Class mutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
:
    public mutURoughWallFunctionFvPatchScalarField
{
public:

    //- Runtime type information
    TypeName("mutSpalartAllmarasStandardRoughWallFunction");


    // Constructors

        //- Construct from patch and internal field
        mutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
        //  onto a new patch
        mutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
        (
            const mutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
        (
            const mutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        mutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
        (
            const mutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField&,
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
                new mutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
                (
                   *this,
                    iF
                )
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
