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
    Foam::MapFvVolField

Description
    Map volume internal field on topology change.  This is a partial
    template specialisation, see MapGeometricFields.

\*---------------------------------------------------------------------------*/

#ifndef MapFvVolField_H
#define MapFvVolField_H

#include "Field.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class MeshMapper>
class MapInternalField<Type, MeshMapper, volMesh>
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


template<class Type, class MeshMapper>
void MapInternalField<Type, MeshMapper, volMesh>::operator()
(
    Field<Type>& field,
    const MeshMapper& mapper
) const
{
    if (field.size() != mapper.volMap().sizeBeforeMapping())
    {
        FatalErrorIn
        (
            "void MapInternalField<Type, MeshMapper, volMesh>::operator()\n"
            "(\n"
            "    Field<Type>& field,\n"
            "    const MeshMapper& mapper\n"
            ") const"
        )  << "Incompatible size before mapping.  Field size: " << field.size()
           << " map size: " << mapper.volMap().sizeBeforeMapping()
           << abort(FatalError);
    }

    field.autoMap(mapper.volMap());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
