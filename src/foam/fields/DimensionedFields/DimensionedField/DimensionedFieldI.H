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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
inline const DimensionedField<Type, GeoMesh>&
DimensionedField<Type, GeoMesh>::null()
{
    // Dangerous, but no null constructor.  Reconsider.  HJ, 10/Apr/2013
    return *reinterpret_cast< DimensionedField<Type, GeoMesh>* >(0);
}


template<class Type, class GeoMesh>
inline const typename GeoMesh::Mesh&
DimensionedField<Type, GeoMesh>::mesh() const
{
    return mesh_;
}


template<class Type, class GeoMesh>
inline const dimensionSet& DimensionedField<Type, GeoMesh>::dimensions() const
{
    return dimensions_;
}

template<class Type, class GeoMesh>
inline dimensionSet& DimensionedField<Type, GeoMesh>::dimensions()
{
    return dimensions_;
}


template<class Type, class GeoMesh>
inline const Field<Type>& DimensionedField<Type, GeoMesh>::field() const
{
    return *this;
}

template<class Type, class GeoMesh>
inline Field<Type>& DimensionedField<Type, GeoMesh>::field()
{
    return *this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
