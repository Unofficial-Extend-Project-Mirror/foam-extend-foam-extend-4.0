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

\*---------------------------------------------------------------------------*/

#include "backwardsCompatibilityWallFunctions.H"

#include "calculatedFvPatchField.H"
#include "alphatWallFunctionFvPatchScalarField.H"
#include "mutWallFunctionFvPatchScalarField.H"
#include "mutLowReWallFunctionFvPatchScalarField.H"
#include "epsilonWallFunctionFvPatchScalarField.H"
#include "kqRWallFunctionFvPatchField.H"
#include "omegaWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> autoCreateAlphat
(
    const word& fieldName,
    const fvMesh& mesh,
    const objectRegistry& obj
)
{
    IOobject alphatHeader
    (
        fieldName,
        mesh.time().timeName(),
        obj,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (alphatHeader.headerOk())
    {
        return tmp<volScalarField>(new volScalarField(alphatHeader, mesh));
    }
    else
    {
        Info<< "--> Creating " << fieldName
            << " to employ run-time selectable wall functions" << endl;

        const fvBoundaryMesh& bm = mesh.boundary();

        wordList alphatBoundaryTypes(bm.size());

        forAll(bm, patchI)
        {
            if (bm[patchI].isWall())
            {
                alphatBoundaryTypes[patchI] =
                    RASModels::alphatWallFunctionFvPatchScalarField::typeName;
            }
            else
            {
                alphatBoundaryTypes[patchI] =
                    calculatedFvPatchField<scalar>::typeName;
            }
        }

        tmp<volScalarField> alphat
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    obj,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar("zero", dimDensity*dimArea/dimTime, 0.0),
                alphatBoundaryTypes
            )
        );

        Info<< "    Writing new " << fieldName << endl;
        alphat().write();

        return alphat;
    }
}


tmp<volScalarField> autoCreateAlphat
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return autoCreateAlphat(fieldName, mesh, mesh);
}


tmp<volScalarField> autoCreateMut
(
    const word& fieldName,
    const fvMesh& mesh,
    const objectRegistry& obj
)
{
    IOobject mutHeader
    (
        fieldName,
        mesh.time().timeName(),
        obj,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (mutHeader.headerOk())
    {
        return tmp<volScalarField>(new volScalarField(mutHeader, mesh));
    }
    else
    {
        Info<< "--> Creating " << fieldName
            << " to employ run-time selectable wall functions" << endl;

        const fvBoundaryMesh& bm = mesh.boundary();

        wordList mutBoundaryTypes(bm.size());

        forAll(bm, patchI)
        {
            if (bm[patchI].isWall())
            {
                mutBoundaryTypes[patchI] =
                    RASModels::mutWallFunctionFvPatchScalarField::typeName;
            }
            else
            {
                mutBoundaryTypes[patchI] =
                    calculatedFvPatchField<scalar>::typeName;
            }
        }

        tmp<volScalarField> mut
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    obj,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar("zero", dimDensity*dimArea/dimTime, 0.0),
                mutBoundaryTypes
            )
        );

        Info<< "    Writing new " << fieldName << endl;
        mut().write();

        return mut;
    }
}


tmp<volScalarField> autoCreateMut
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return autoCreateMut(fieldName, mesh, mesh);
}


tmp<volScalarField> autoCreateLowReMut
(
    const word& fieldName,
    const fvMesh& mesh,
    const objectRegistry& obj
)
{
    IOobject mutHeader
    (
        fieldName,
        mesh.time().timeName(),
        obj,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (mutHeader.headerOk())
    {
        return tmp<volScalarField>(new volScalarField(mutHeader, mesh));
    }
    else
    {
        Info<< "--> Creating " << fieldName
            << " to employ run-time selectable wall functions" << endl;

        const fvBoundaryMesh& bm = mesh.boundary();

        wordList mutBoundaryTypes(bm.size());

        forAll(bm, patchI)
        {
            if (bm[patchI].isWall())
            {
                mutBoundaryTypes[patchI] =
                    RASModels::mutLowReWallFunctionFvPatchScalarField::typeName;
            }
            else
            {
                mutBoundaryTypes[patchI] =
                    calculatedFvPatchField<scalar>::typeName;
            }
        }

        tmp<volScalarField> mut
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    obj,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar("zero", dimDensity*dimArea/dimTime, 0.0),
                mutBoundaryTypes
            )
        );

        Info<< "    Writing new " << fieldName << endl;
        mut().write();

        return mut;
    }
}


tmp<volScalarField> autoCreateLowReMut
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return autoCreateLowReMut(fieldName, mesh, mesh);
}


tmp<volScalarField> autoCreateEpsilon
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
            RASModels::epsilonWallFunctionFvPatchScalarField
        >
        (
            fieldName,
            mesh,
            obj
        );
}


tmp<volScalarField> autoCreateEpsilon
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            RASModels::epsilonWallFunctionFvPatchScalarField
        >
        (
            fieldName,
            mesh,
            mesh
        );
}


tmp<volScalarField> autoCreateOmega
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
            RASModels::omegaWallFunctionFvPatchScalarField
        >
        (
            fieldName,
            mesh,
            obj
        );
}


tmp<volScalarField> autoCreateOmega
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            RASModels::omegaWallFunctionFvPatchScalarField
        >
        (
            fieldName,
            mesh,
            mesh
        );
}


tmp<volScalarField> autoCreateK
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
            RASModels::kqRWallFunctionFvPatchField<scalar>
        >
        (
            fieldName,
            mesh,
            obj
        );
}


tmp<volScalarField> autoCreateK
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            RASModels::kqRWallFunctionFvPatchField<scalar>
        >
        (
            fieldName,
            mesh,
            mesh
        );
}


tmp<volScalarField> autoCreateQ
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
            RASModels::kqRWallFunctionFvPatchField<scalar>
        >
        (
            fieldName,
            mesh,
            obj
        );
}


tmp<volScalarField> autoCreateQ
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            RASModels::kqRWallFunctionFvPatchField<scalar>
        >
        (
            fieldName,
            mesh,
            mesh
        );
}


tmp<volSymmTensorField> autoCreateR
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
            RASModels::kqRWallFunctionFvPatchField<symmTensor>
        >
        (
            fieldName,
            mesh,
            obj
        );
}


tmp<volSymmTensorField> autoCreateR
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            symmTensor,
            RASModels::kqRWallFunctionFvPatchField<symmTensor>
        >
        (
            fieldName,
            mesh,
            mesh
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //

