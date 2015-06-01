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

#include "sampledTriSurfaceMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledTriSurfaceMesh::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    // One value per face
    tmp<Field<Type> > tvalues(new Field<Type>(cellLabels_.size()));
    Field<Type>& values = tvalues();

    forAll(cellLabels_, triI)
    {
        values[triI] = vField[cellLabels_[triI]];
    }

    return tvalues;
}


template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledTriSurfaceMesh::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    // One value per vertex
    tmp<Field<Type> > tvalues(new Field<Type>(pointToFace_.size()));
    Field<Type>& values = tvalues();

    forAll(pointToFace_, pointI)
    {
        label triI = pointToFace_[pointI];
        label cellI = cellLabels_[triI];

        values[pointI] = interpolator.interpolate(points()[pointI], cellI);
    }

    return tvalues;
}


// ************************************************************************* //
