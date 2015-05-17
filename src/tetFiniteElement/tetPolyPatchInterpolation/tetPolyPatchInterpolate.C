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

Description

\*---------------------------------------------------------------------------*/

#include "tetPolyPatchInterpolation.H"
#include "faceTetPolyPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::tetPolyPatchInterpolation::faceToPointInterpolate
(
    const Field<Type>& ff
) const
{
    if (ff.size() != patch_.patch().size())
    {
        FatalErrorIn
        (
            "tmp<Foam::Field<Type> >\n"
            "tetPolyPatchInterpolation::faceToPointInterpolate\n"
            "(\n"
            "    const Field<Type>& ff\n"
            ") const"
        )  << "Field size: " << ff.size() << " does not match number of faces: "
           << patch_.patch().size()
           << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>(patch_.size())
    );
    Field<Type>& result = tresult();

    // Insert the point values first
    label i = 0;

    const Field<Type> pointRes = interpolator_.faceToPointInterpolate(ff);

    forAll (pointRes, pointI)
    {
        result[i] = pointRes[pointI];
        i++;
    }

    // Insert the face centre values; no interpolation necessary
    forAll (ff, faceI)
    {
        result[i] = ff[faceI];
        i++;
    }

    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::tetPolyPatchInterpolation::faceToPointInterpolate
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = this->faceToPointInterpolate(tff());
    tff.clear();
    return tint;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::tetPolyPatchInterpolation::pointToPointInterpolate
(
    const Field<Type>& ff
) const
{
    if (ff.size() != patch_.patch().nPoints())
    {
        FatalErrorIn
        (
            "tmp<Foam::Field<Type> >\n"
            "tetPolyPatchInterpolation::pointToPointInterpolate\n"
            "(\n"
            "    const Field<Type>& ff\n"
            ") const"
        )  << "Field size: " << ff.size()
           << " does not match number of points: "
           << patch_.patch().nPoints()
           << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>(patch_.size())
    );
    Field<Type>& result = tresult();

    // Insert the point values first; no interpolation necessary
    label i = 0;

    forAll (ff, pointI)
    {
        result[i] = ff[pointI];
        i++;
    }

    // Insert the face centre values
    const Field<Type> faceRes = interpolator_.pointToFaceInterpolate(ff);

    forAll (faceRes, faceI)
    {
        result[i] = faceRes[faceI];
        i++;
    }

    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::tetPolyPatchInterpolation::pointToPointInterpolate
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = this->pointToPointInterpolate(tff());
    tff.clear();
    return tint;
}


// ************************************************************************* //
