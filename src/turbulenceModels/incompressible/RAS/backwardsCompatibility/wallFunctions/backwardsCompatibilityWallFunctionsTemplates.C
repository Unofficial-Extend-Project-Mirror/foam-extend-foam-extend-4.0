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
#include "foamTime.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class PatchType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
autoCreateWallFunctionField
(
    const word& fieldName,
    const fvMesh& mesh,
    const objectRegistry& obj
)
{
    IOobject nutHeader
    (
        "nut",
        mesh.time().timeName(),
        obj,
        IOobject::MUST_READ
    );

    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (nutHeader.headerOk())
    {
        return tmp<fieldType>
        (
            new fieldType
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    obj,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh
            )
        );
    }
    else
    {
        Info<< "--> Upgrading " << fieldName
            << " to employ run-time selectable wall functions" << endl;

        // Read existing field
        IOobject ioObj
        (
            fieldName,
            mesh.time().timeName(),
            obj,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        tmp<fieldType> fieldOrig
        (
            new fieldType
            (
                ioObj,
                mesh
            )
        );

        // rename file
        Info<< "    Backup original " << fieldName << " to "
            << fieldName << ".old" << endl;
        mvBak(ioObj.objectPath(), "old");


        PtrList<fvPatchField<Type> > newPatchFields(mesh.boundary().size());

        forAll(newPatchFields, patchI)
        {
            if (mesh.boundary()[patchI].isWall())
            {
                newPatchFields.set
                (
                    patchI,
                    new PatchType
                    (
                        mesh.boundary()[patchI],
                        fieldOrig().dimensionedInternalField()
                    )
                );
                newPatchFields[patchI] == fieldOrig().boundaryField()[patchI];
            }
            else
            {
                newPatchFields.set
                (
                    patchI,
                    fieldOrig().boundaryField()[patchI].clone()
                );
            }
        }

        tmp<fieldType> fieldNew
        (
            new fieldType
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
                fieldOrig().dimensions(),
                fieldOrig().internalField(),
                newPatchFields
            )
        );

        Info<< "    Writing updated " << fieldName << endl;
        fieldNew().write();

        return fieldNew;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
