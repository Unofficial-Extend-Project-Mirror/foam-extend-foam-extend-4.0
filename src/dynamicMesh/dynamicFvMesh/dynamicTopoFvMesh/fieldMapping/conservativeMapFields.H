/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "topoMapper.H"
#include "topoCellMapper.H"
#include "topoSurfaceMapper.H"
#include "topoBoundaryMeshMapper.H"

namespace Foam
{

// Conservatively map all volFields in the registry
template <class Type>
void conservativeMapVolFields(const topoMapper& mapper)
{
    // Define a few typedefs for convenience
    typedef typename outerProduct<vector, Type>::type gCmptType;
    typedef GeometricField<Type, fvPatchField, volMesh> volType;
    typedef GeometricField<gCmptType, fvPatchField, volMesh> gradVolType;

    HashTable<const volType*> fields(mapper.mesh().lookupClass<volType>());

    // Store old-times before mapping
    forAllIter(typename HashTable<const volType*>, fields, fIter)
    {
        fIter()->storeOldTimes();
    }

    // Fetch internal/boundary mappers
    const topoCellMapper& fMap = mapper.volMap();
    const topoBoundaryMeshMapper& bMap = mapper.boundaryMap();

    // Now map all fields
    forAllIter(typename HashTable<const volType*>, fields, fIter)
    {
        volType& field = const_cast<volType&>(*fIter());

        if (fvMesh::debug)
        {
            Info<< "Conservatively mapping "
                << field.typeName
                << ' ' << field.name()
                << endl;
        }

        // Map the internal field
        fMap.mapInternalField
        (
            field.name(),
            mapper.gradient<gradVolType>(field.name()).internalField(),
            field.internalField()
        );

        // Map patch fields
        forAll(bMap, patchI)
        {
            bMap[patchI].mapPatchField
            (
                field.name(),
                field.boundaryField()[patchI]
            );
        }

        // Set the field instance
        field.instance() = field.mesh().thisDb().time().timeName();
    }
}


// Conservatively map all surfaceFields in the registry
template <class Type>
void conservativeMapSurfaceFields(const topoMapper& mapper)
{
    // Define a few typedefs for convenience
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> surfType;

    HashTable<const surfType*> fields(mapper.mesh().lookupClass<surfType>());

    // Store old-times before mapping
    forAllIter(typename HashTable<const surfType*>, fields, fIter)
    {
        fIter()->storeOldTimes();
    }

    // Fetch internal/boundary mappers
    const topoSurfaceMapper& fMap = mapper.surfaceMap();
    const topoBoundaryMeshMapper& bMap = mapper.boundaryMap();

    // Now map all fields
    forAllIter(typename HashTable<const surfType*>, fields, fIter)
    {
        surfType& field = const_cast<surfType&>(*fIter());

        if (fvMesh::debug)
        {
            Info<< "Conservatively mapping "
                << field.typeName
                << ' ' << field.name()
                << endl;
        }

        // Map the internal field
        fMap.mapInternalField
        (
            field.name(),
            field.internalField()
        );

        // Map patch fields
        forAll(bMap, patchI)
        {
            bMap[patchI].mapPatchField
            (
                field.name(),
                field.boundaryField()[patchI]
            );
        }

        // Set the field instance
        field.instance() = field.mesh().thisDb().time().timeName();
    }
}


} // End namespace Foam

// ************************************************************************* //
