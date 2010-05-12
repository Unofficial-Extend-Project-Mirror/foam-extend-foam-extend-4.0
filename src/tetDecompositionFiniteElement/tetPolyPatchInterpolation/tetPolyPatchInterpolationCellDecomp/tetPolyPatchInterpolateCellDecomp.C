/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/

#include "tetPolyPatchInterpolationCellDecomp.H"
#include "faceTetPolyPatchCellDecomp.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::tetPolyPatchInterpolationCellDecomp::faceToPointInterpolate
(
    const Field<Type>& ff
) const
{
    if (ff.size() != patch_.patch().size())
    {
        FatalErrorIn
        (
            "tmp<Foam::Field<Type> >\n"
            "tetPolyPatchInterpolationCellDecomp::faceToPointInterpolate\n"
            "(\n"
            "    const Field<Type>& ff\n"
            ") const"
        )  << "Field size: " << ff.size() << " does not match number of faces: "
           << patch_.patch().size()
           << abort(FatalError);
    }

    return interpolator_.faceToPointInterpolate(ff);
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::tetPolyPatchInterpolationCellDecomp::faceToPointInterpolate
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
Foam::tetPolyPatchInterpolationCellDecomp::pointToPointInterpolate
(
    const Field<Type>& ff
) const
{
    if (ff.size() != patch_.patch().nPoints())
    {
        FatalErrorIn
        (
            "tmp<Foam::Field<Type> >\n"
            "tetPolyPatchInterpolationCellDecomp::pointToPointInterpolate\n"
            "(\n"
            "    const Field<Type>& ff\n"
            ") const"
        )  << "Field size: " << ff.size()
           << " does not match number of points: "
           << patch_.patch().nPoints()
           << abort(FatalError);
    }

    return tmp<Field<Type> >(ff);
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::tetPolyPatchInterpolationCellDecomp::pointToPointInterpolate
(
    const tmp<Field<Type> >& tff
) const
{
  if (tff().size() != patch_.patch().nPoints())
    {
        FatalErrorIn
        (
            "tmp<Foam::Field<Type> >\n"
            "tetPolyPatchInterpolationCellDecomp::pointToPointInterpolate\n"
            "(\n"
            "    const tmp<Field<Type> >\n"
            ") const"
        )  << "Field size: " << tff().size()
           << " does not match number of points: "
           << patch_.patch().nPoints()
           << abort(FatalError);
    }

    return tff;
}


// ************************************************************************* //
