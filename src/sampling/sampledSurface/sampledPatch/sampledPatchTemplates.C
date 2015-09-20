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

#include "sampledPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledPatch::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    // One value per face
    tmp<Field<Type> > tvalues(new Field<Type>(patchFaceLabels_.size()));
    Field<Type>& values = tvalues();

    if (patchIndex() != -1)
    {
        const Field<Type>& bField = vField.boundaryField()[patchIndex()];

        forAll(patchFaceLabels_, elemI)
        {
            values[elemI] = bField[patchFaceLabels_[elemI]];
        }
    }

    return tvalues;
}


template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledPatch::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    // One value per vertex
    tmp<Field<Type> > tvalues(new Field<Type>(points().size()));
    Field<Type>& values = tvalues();

    if (patchIndex() != -1)
    {
        const polyPatch& patch = mesh().boundaryMesh()[patchIndex()];
        const labelList& own = mesh().faceOwner();

        boolList pointDone(points().size(), false);

        forAll(faces(), cutFaceI)
        {
            const face& f = faces()[cutFaceI];

            forAll(f, faceVertI)
            {
                label pointI = f[faceVertI];

                if (!pointDone[pointI])
                {
                    label faceI = patchFaceLabels()[cutFaceI] + patch.start();
                    label cellI = own[faceI];

                    values[pointI] = interpolator.interpolate
                    (
                        points()[pointI],
                        cellI,
                        faceI
                    );
                    pointDone[pointI] = true;
                }
            }
        }
    }

    return tvalues;
}


// ************************************************************************* //
