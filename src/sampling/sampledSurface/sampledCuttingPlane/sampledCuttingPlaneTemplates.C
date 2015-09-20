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

#include "volPointInterpolation.H"
#include "sampledCuttingPlane.H"
#include "isoSurface.H"
#include "volFieldsFwd.H"
#include "pointFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledCuttingPlane::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    return tmp<Field<Type> >(new Field<Type>(vField, surface().meshCells()));
}


template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledCuttingPlane::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    // Get fields to sample. Assume volPointInterpolation!
    const GeometricField<Type, fvPatchField, volMesh>& volFld =
        interpolator.psi();

    if (subMeshPtr_.valid())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh> > tvolSubFld =
            subMeshPtr_().interpolate(volFld);

        const GeometricField<Type, fvPatchField, volMesh>& volSubFld =
            tvolSubFld();

        tmp<GeometricField<Type, pointPatchField, pointMesh> > tpointSubFld =
            volPointInterpolation::New(volSubFld.mesh()).interpolate(volSubFld);

        // Sample.
        return surface().interpolate(volSubFld, tpointSubFld());
    }
    else
    {
        tmp<GeometricField<Type, pointPatchField, pointMesh> > tpointFld
        (
            volPointInterpolation::New(volFld.mesh()).interpolate(volFld)
        );

        // Sample.
        return surface().interpolate(volFld, tpointFld());
    }
}


// ************************************************************************* //
