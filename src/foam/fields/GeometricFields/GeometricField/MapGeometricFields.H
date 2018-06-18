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
    Foam::MapInternalField
    Foam::MapGeometricFields

Description
    Generic internal field mapper.  For "real" mapping, add template
    specialisations for mapping of internal fields depending on mesh
    type.

\*---------------------------------------------------------------------------*/

#ifndef MapGeometricFields_H
#define MapGeometricFields_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class MeshMapper, class GeoMesh>
class MapInternalField
{
public:

    MapInternalField()
    {}

    void operator()
    (
        Field<Type>& field,
        const MeshMapper& mapper
    ) const;
};


//- Generic Geometric field mapper.
//  For "real" mapping, add template specialisations
//  for mapping of internal fields depending on mesh type.
template
<
    class Type,
    template<class> class PatchField,
    class MeshMapper,
    class GeoMesh
>
void MapGeometricFields
(
    const MeshMapper& mapper
)
{
    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        mapper.thisDb().objectRegistry::template lookupClass
            <GeometricField<Type, PatchField, GeoMesh> >()
    );

    // It is necessary to enforce that all old-time fields are stored
    // before the mapping is performed.  Otherwise, if the
    // old-time-level field is mapped before the field itself, sizes
    // will not match.

    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, PatchField, GeoMesh>& field =
            const_cast<GeometricField<Type, PatchField, GeoMesh>&>
            (*fieldIter());

        // Note: check can be removed once pointFields are actually stored on
        //       the pointMesh instead of now on the polyMesh!
        if (&field.mesh() == &mapper.mesh())
        {
            field.storeOldTimes();
        }
    }

    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, PatchField, GeoMesh>& field =
            const_cast<GeometricField<Type, PatchField, GeoMesh>&>
            (*fieldIter());

        if (&field.mesh() == &mapper.mesh())
        {
            if (polyMesh::debug)
            {
                InfoIn("void MapGeometricFields(const MeshMapper& mapper)")
                    << "Mapping " << field.typeName << ' ' << field.name()
                    << endl;
            }

            // Map the internal field
            MapInternalField<Type, MeshMapper, GeoMesh>()
            (
                field.internalField(),
                mapper
            );

            // Repatch boundary if needed
            if
            (
                field.boundaryField().size()
             != mapper.mesh().boundary().size()
            )
            {
                if (polyMesh::debug)
                {
                    InfoIn("void MapGeometricFields(const MeshMapper& mapper)")
                        << "Resizing boundary field for "
                        << field.typeName << ' ' << field.name()
                        << endl;
                }

                field.boundaryField().setSize(mapper.mesh().boundary().size());
            }

            // Get repatch flag
            const boolList& resetPatchFlag = mapper.resetPatchFlag();

            if (resetPatchFlag.size() != field.boundaryField().size())
            {
                FatalErrorIn
                (
                    "void MapGeometricFields(const MeshMapper& mapper)"
                )   << "Incorrect resetPatchFlag array size.  Boundary: "
                    << field.boundaryField().size() << " flag: "
                    << resetPatchFlag.size()
                    << abort(FatalError);
            }

            // Map the patch fields
            forAll (field.boundaryField(), patchI)
            {
                // Flag reset can only take place for the constrained patches
                // such as empty or processor
                if (resetPatchFlag[patchI])
                {
                    // Build a new patchField if reset is true
                    field.boundaryField().set
                    (
                        patchI,
                        PatchField<Type>::New
                        (
                            mapper.mesh().boundary()[patchI].type(),
                            field.mesh().boundary()[patchI],
                            field.dimensionedInternalField()
                        )
                    );

                    // Set to zero to avoid NaNs?
                }
                else
                {
                    // No reset: auto-map

                    // Cannot check sizes for patch fields because of
                    // empty fields in FV and because point fields get
                    // their size from the patch which has already been resized

                    field.boundaryField()[patchI].autoMap
                    (
                        mapper.boundaryMap()[patchI]
                    );
                }
            }

            field.instance() = field.time().timeName();
        }
        else
        {
            if (polyMesh::debug)
            {
                InfoIn("void MapGeometricFields(const MeshMapper& mapper)")
                    << "Not mapping " << field.typeName << ' ' << field.name()
                    << " since originating mesh differs from that of mapper."
                    << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
