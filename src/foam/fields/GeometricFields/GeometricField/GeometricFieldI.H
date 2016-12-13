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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
inline const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::null()
{
    return *reinterpret_cast< GeometricField<Type, PatchField, GeoMesh>* >(0);
}


template<class Type, template<class> class PatchField, class GeoMesh>
inline
const typename
Foam::GeometricField<Type, PatchField, GeoMesh>::DimensionedInternalField&
Foam::GeometricField<Type, PatchField, GeoMesh>::
dimensionedInternalField() const
{
    return *this;
}


template<class Type, template<class> class PatchField, class GeoMesh>
inline
const typename
Foam::GeometricField<Type, PatchField, GeoMesh>::InternalField&
Foam::GeometricField<Type, PatchField, GeoMesh>::internalField() const
{
    return *this;
}


template<class Type, template<class> class PatchField, class GeoMesh>
inline const typename Foam::GeometricField<Type, PatchField, GeoMesh>::
GeometricBoundaryField&
Foam::GeometricField<Type, PatchField, GeoMesh>::boundaryField() const
{
    return boundaryField_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
inline Foam::label
Foam::GeometricField<Type, PatchField, GeoMesh>::timeIndex() const
{
    return timeIndex_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
inline Foam::label&
Foam::GeometricField<Type, PatchField, GeoMesh>::timeIndex()
{
    return timeIndex_;
}


// ************************************************************************* //
