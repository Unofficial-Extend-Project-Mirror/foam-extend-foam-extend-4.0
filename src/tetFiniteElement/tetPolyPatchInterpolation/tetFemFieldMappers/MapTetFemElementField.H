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
    MapTetFemElementField

Description
    Map element internal field on topology change.  This is a partial
    template specialisation, see MapGeometricFields.

\*---------------------------------------------------------------------------*/

#ifndef MapTetFemElementField_H
#define MapTetFemElementField_H

#include "Field.H"
#include "elementMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class MeshMapper>
class MapInternalField<Type, MeshMapper, elementMesh>
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
void MapInternalField<Type, MeshMapper, elementMesh>::operator()
(
    Field<Type>& field,
    const MeshMapper& mapper
) const
{
    if (field.size() != mapper.elementMap().sizeBeforeMapping())
    {
        FatalErrorIn
        (
            "void MapInternalField<Type, MeshMapper, elementMesh>::operator()\n"
            "(\n"
            "    Field<Type>& field,\n"
            "    const MeshMapper& mapper\n"
            ") const"
        )  << "Incompatible size before mapping.  Field size: " << field.size()
           << " map size: " << mapper.elementMap().sizeBeforeMapping()
           << abort(FatalError);
    }

    field.autoMap(mapper.elementMap());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
