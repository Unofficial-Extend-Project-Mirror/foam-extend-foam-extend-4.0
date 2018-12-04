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

\*---------------------------------------------------------------------------*/

#include "backwardsCompatibilityIbWallFunctions.H"

#include "calculatedFvPatchField.H"
#include "nutkWallFunctionFvPatchScalarField.H"
#include "nutLowReWallFunctionFvPatchScalarField.H"
#include "epsilonWallFunctionFvPatchScalarField.H"
#include "kqRWallFunctionFvPatchField.H"
#include "RWallFunctionFvPatchSymmTensorField.H"
#include "omegaWallFunctionFvPatchScalarField.H"
#include "immersedBoundaryEpsilonWallFunctionFvPatchScalarField.H"
#include "immersedBoundaryKqRWallFunctionFvPatchField.H"
#include "immersedBoundaryOmegaWallFunctionFvPatchScalarField.H"
#include "immersedBoundaryNutWallFunctionFvPatchScalarField.H"
#include "immersedBoundaryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> autoCreateIbNut
(
    const word& fieldName,
    const fvMesh& mesh,
    const objectRegistry& obj
)
{
    IOobject nutHeader
    (
        fieldName,
        mesh.time().timeName(),
        obj,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (nutHeader.headerOk())
    {
        return tmp<volScalarField>(new volScalarField(nutHeader, mesh));
    }
    else
    {
        Info<< "--> Creating " << fieldName
            << " to employ run-time selectable wall functions" << endl;

        const fvBoundaryMesh& bm = mesh.boundary();

        wordList nutBoundaryTypes(bm.size());

        forAll(bm, patchI)
        {
            if (isA<immersedBoundaryFvPatch>(bm[patchI]) && bm[patchI].isWall())
            {
                nutBoundaryTypes[patchI] =
                    RASModels::immersedBoundaryNutWallFunctionFvPatchScalarField::typeName;
            }
            else if (bm[patchI].isWall())
            {
                nutBoundaryTypes[patchI] =
                    RASModels::nutkWallFunctionFvPatchScalarField::typeName;
            }
            else
            {
                nutBoundaryTypes[patchI] =
                    calculatedFvPatchField<scalar>::typeName;
            }
        }

        tmp<volScalarField> nut
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar("zero", dimArea/dimTime, 0.0),
                nutBoundaryTypes
            )
        );

        Info<< "    Writing new " << fieldName << endl;
        nut().write();

        return nut;
    }
}


tmp<volScalarField> autoCreateIbNut
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return autoCreateIbNut(fieldName, mesh, mesh);
}


tmp<volScalarField> autoCreateIbEpsilon
(
    const word& fieldName,
    const fvMesh& mesh,
    const objectRegistry& obj
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            RASModels::epsilonWallFunctionFvPatchScalarField,
            RASModels::immersedBoundaryEpsilonWallFunctionFvPatchScalarField
        >
        (
            fieldName,
            mesh,
            obj
        );
}


tmp<volScalarField> autoCreateIbEpsilon
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            RASModels::epsilonWallFunctionFvPatchScalarField,
            RASModels::immersedBoundaryEpsilonWallFunctionFvPatchScalarField
        >
        (
            fieldName,
            mesh,
            mesh
        );
}


tmp<volScalarField> autoCreateIbOmega
(
    const word& fieldName,
    const fvMesh& mesh,
    const objectRegistry& obj
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            RASModels::omegaWallFunctionFvPatchScalarField,
            RASModels::immersedBoundaryOmegaWallFunctionFvPatchScalarField
        >
        (
            fieldName,
            mesh,
            obj
        );
}


tmp<volScalarField> autoCreateIbOmega
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            RASModels::omegaWallFunctionFvPatchScalarField,
            RASModels::immersedBoundaryOmegaWallFunctionFvPatchScalarField
        >
        (
            fieldName,
            mesh,
            mesh
        );
}


tmp<volScalarField> autoCreateIbK
(
    const word& fieldName,
    const fvMesh& mesh,
    const objectRegistry& obj
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            RASModels::kqRWallFunctionFvPatchField<scalar>,
            RASModels::immersedBoundaryKqRWallFunctionFvPatchField<scalar>
        >
        (
            fieldName,
            mesh,
            obj
        );
}


tmp<volScalarField> autoCreateIbK
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            RASModels::kqRWallFunctionFvPatchField<scalar>,
            RASModels::immersedBoundaryKqRWallFunctionFvPatchField<scalar>
        >
        (
            fieldName,
            mesh,
            mesh
        );
}


tmp<volScalarField> autoCreateIbQ
(
    const word& fieldName,
    const fvMesh& mesh,
    const objectRegistry& obj
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            RASModels::kqRWallFunctionFvPatchField<scalar>,
            RASModels::immersedBoundaryKqRWallFunctionFvPatchField<scalar>
        >
        (
            fieldName,
            mesh,
            obj
        );
}


tmp<volScalarField> autoCreateIbQ
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            RASModels::kqRWallFunctionFvPatchField<scalar>,
            RASModels::immersedBoundaryKqRWallFunctionFvPatchField<scalar>
        >
        (
            fieldName,
            mesh,
            mesh
        );
}


tmp<volSymmTensorField> autoCreateIbR
(
    const word& fieldName,
    const fvMesh& mesh,
    const objectRegistry& obj
)
{
    return
        autoCreateWallFunctionField
        <
            symmTensor,
            RASModels::kqRWallFunctionFvPatchField<symmTensor>,
            RASModels::immersedBoundaryKqRWallFunctionFvPatchField<symmTensor>
        >
        (
            fieldName,
            mesh,
            obj
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

