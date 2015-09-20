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
    Foam::MapFaEdgeField

Description
    Map edge internal field on topology change.  This is a partial
    template specialisation, see MapGeometricFields.

\*---------------------------------------------------------------------------*/

#ifndef MapFaEdgeField_H
#define MapFaEdgeField_H

#include "Field.H"
#include "edgeMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class MeshMapper>
class MapInternalField<Type, MeshMapper, edgeMesh>
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
void MapInternalField<Type, MeshMapper, edgeMesh>::operator()
(
    Field<Type>& field,
    const MeshMapper& mapper
) const
{
    if (field.size() != mapper.edgeMap().sizeBeforeMapping())
    {
        FatalErrorIn
        (
            "void MapInternalField<Type, MeshMapper, edgeMesh>::operator()\n"
            "(\n"
            "    Field<Type>& field,\n"
            "    const MeshMapper& mapper\n"
            ") const"
        )  << "Incompatible size before mapping.  Field size: " << field.size()
           << " map size: " << mapper.edgeMap().sizeBeforeMapping()
           << abort(FatalError);
    }

    Info<< "Note: No mapping of FA edge fields.  Please reconsider algorithm"
        << endl;

    field.autoMap(mapper.edgeMap());

    // Flip the flux
//     const labelList flipFaces = mapper.edgeMap().flipFaceFlux().toc();

//     forAll (flipFaces, i)
//     {
//         if (flipFaces[i] < field.size())
//         {
//             field[flipFaces[i]] *= -1.0;
//         }
//     }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
