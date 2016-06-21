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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Fethi Tekin, All rights reserved.
    Oliver Borm, All rights reserved.

\*---------------------------------------------------------------------------*/

#include "overlapGgiPolyPatch.H"
#include "RodriguesRotation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::overlapGgiPolyPatch::expandData(const Field<Type>& pf) const
{
    // Check and expand the field from patch size to zone size
    if (pf.size() != size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > overlapGgiPolyPatch::expandData"
            "("
            "    const Field<Type>& pf"
            ") const"
        )   << "Incorrect patch field size.  Field size: "
            << pf.size() << " patch size: " << size()
            << abort(FatalError);
    }

    const label ncp = nCopies();

    const scalar myAngle = 360.0/scalar(ncp);

    tmp<Field<Type> > texpandField
    (
        new Field<Type>(ncp*zone().size(), pTraits<Type>::zero)
    );

    Field<Type>& expandField = texpandField();

    for (label copyI = 0; copyI < ncp; copyI++)
    {
        // Calculate transform
        const tensor curRotation =
            RodriguesRotation(rotationAxis_, copyI*myAngle);

        const label offset = copyI*zone().size();

        forAll (pf, faceI)
        {
            const label zId = zone().whichFace(start() + faceI);
            expandField[offset + zId] = transform(curRotation, pf[faceI]);
        }
    }

    if (!localParallel())
    {
        reduce(expandField, sumOp<Field<Type> >());
    }

    return texpandField;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::overlapGgiPolyPatch::interpolate(const Field<Type>& pf) const
{
    // Check and expand the field from patch size to zone size
    if (pf.size() != shadow().size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > ggiPolyPatch::interpolate"
            "("
            "    const Field<Type>& pf"
            ") const"
        )   << "Incorrect slave patch field size.  Field size: "
            << pf.size() << " patch size: " << shadow().size()
            << abort(FatalError);
    }

    // Expand data
    tmp<Field<Type> > expanddata = shadow().expandData(pf);

    tmp<Field<Type> > tresult(new Field<Type>());
    Field<Type>& result = tresult();

    if (master())
    {
        // Expand slave data
        result = patchToPatch().slaveToMaster(expanddata);
    }
    else
    {
        // Expand master data
        result = patchToPatch().masterToSlave(expanddata);
    }

    // Truncate to size
    result.setSize(size());

    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::overlapGgiPolyPatch::interpolate(const tmp<Field<Type> >& tpf) const
{
    tmp<Field<Type> > tint = interpolate(tpf());
    tpf.clear();
    return tint;
}


// ************************************************************************* //
