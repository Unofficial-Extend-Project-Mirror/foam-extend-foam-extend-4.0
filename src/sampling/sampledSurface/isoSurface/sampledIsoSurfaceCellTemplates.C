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

#include "sampledIsoSurfaceCell.H"
#include "isoSurface.H"
#include "volFieldsFwd.H"
#include "pointFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledIsoSurfaceCell::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    // Recreate geometry if time has changed
    updateGeometry();

    return tmp<Field<Type> >(new Field<Type>(vField, meshCells_));
}


template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledIsoSurfaceCell::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    // Recreate geometry if time has changed
    updateGeometry();

    // One value per point
    tmp<Field<Type> > tvalues(new Field<Type>(points().size()));
    Field<Type>& values = tvalues();

    boolList pointDone(points().size(), false);

    forAll(faces(), cutFaceI)
    {
        const face& f = faces()[cutFaceI];

        forAll(f, faceVertI)
        {
            label pointI = f[faceVertI];

            if (!pointDone[pointI])
            {
                values[pointI] = interpolator.interpolate
                (
                    points()[pointI],
                    meshCells_[cutFaceI]
                );
                pointDone[pointI] = true;
            }
        }
    }

    return tvalues;
}


// ************************************************************************* //
