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
    Foam::incompressible

Description
    Auto creation of fields to provide backwards compatibility with
    runtime selectable wall functions

SourceFiles
    backwardsCompatibilityIbWallFunctions.C
    backwardsCompatibilityIbWallFunctionsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef backwardsCompatibilityIbWallFunctions_H
#define backwardsCompatibilityIbWallFunctions_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
    //- nut
    tmp<volScalarField> autoCreateIbNut
    (
        const word& fieldName,
        const fvMesh& mesh,
        const objectRegistry& obj
    );

    tmp<volScalarField> autoCreateIbNut
    (
        const word& fieldName,
        const fvMesh& mesh
    );

    //- epsilon
    tmp<volScalarField> autoCreateIbEpsilon
    (
        const word& fieldName,
        const fvMesh& mesh,
        const objectRegistry& obj
    );

    tmp<volScalarField> autoCreateIbEpsilon
    (
        const word& fieldName,
        const fvMesh& mesh
    );

    //- omega
    tmp<volScalarField> autoCreateIbOmega
    (
        const word& fieldName,
        const fvMesh& mesh,
        const objectRegistry& obj
    );

    tmp<volScalarField> autoCreateIbOmega
    (
        const word& fieldName,
        const fvMesh& mesh
    );

    //- k
    tmp<volScalarField> autoCreateIbK
    (
        const word& fieldName,
        const fvMesh& mesh,
        const objectRegistry& obj
    );

    tmp<volScalarField> autoCreateIbK
    (
        const word& fieldName,
        const fvMesh& mesh
    );

    //- Q
    tmp<volScalarField> autoCreateIbQ
    (
        const word& fieldName,
        const fvMesh& mesh,
        const objectRegistry& obj
    );

    tmp<volScalarField> autoCreateIbQ
    (
        const word& fieldName,
        const fvMesh& mesh
    );

    //- R
    tmp<volSymmTensorField> autoCreateIbR
    (
        const word& fieldName,
        const fvMesh& mesh,
        const objectRegistry& obj
    );

    tmp<volSymmTensorField> autoCreateIbR
    (
        const word& fieldName,
        const fvMesh& mesh
    );

    //- Helper function to create the new field
    template<class Type, class PatchType, class ibPatchName >
    tmp<GeometricField<Type, fvPatchField, volMesh> >
    autoCreateIbWallFunctionField
    (
        const word& fieldName,
        const fvMesh& mesh,
        const objectRegistry& obj
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "backwardsCompatibilityIbWallFunctionsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
