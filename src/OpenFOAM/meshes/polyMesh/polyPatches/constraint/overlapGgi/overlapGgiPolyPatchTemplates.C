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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Fethi Tekin, All rights reserved.

\*---------------------------------------------------------------------------*/

#include "overlapGgiPolyPatch.H"
#include "RodriguesRotation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::overlapGgiPolyPatch::expandSlaveData(const Field<Type>& spf) const
{
    const scalar slaveAngle = shadow().angle();

    const label ncp = shadow().nCopies();

    tmp<Field<Type> > tef(new Field<Type>(ncp*spf.size()));
    Field<Type>& ef = tef();

    label nFaces = 0;

    for (label copyI = 0; copyI < ncp; copyI++)
    {
        // Calculate transform
        const tensor curRotation =
            RodriguesRotation(rotationAxis_,  copyI*slaveAngle);

        forAll (spf, faceI)
        {
            ef[nFaces] = transform(curRotation, spf[faceI]);
            nFaces++;
        }
    }

    return tef;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::overlapGgiPolyPatch::expandMasterData(const Field<Type>& spf) const
{
    const scalar masterAngle = angle();

    const label ncpm = nCopies();

    tmp<Field<Type> > tef(new Field<Type>(ncpm*spf.size()));
    Field<Type>& ef = tef();

    label nFaces = 0;

    for (label copyI = 0; copyI < ncpm; copyI++)
    {
        // Calculate transform
        const tensor curRotation =
            RodriguesRotation(rotationAxis_, copyI*(masterAngle));

        forAll (spf, faceI)
        {
            ef[nFaces] = transform(curRotation, spf[faceI]);
            nFaces++;
        }
    }

    return tef;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::overlapGgiPolyPatch::interpolate(const Field<Type>& pf) const
{
    if (master())
    {
        // Expand slave data
        tmp<Field<Type> > expandslave = expandSlaveData(pf);

        tmp<Field<Type> > tresult= patchToPatch().slaveToMaster(expandslave);
        // Truncate to size
        tresult().setSize(size());

        return tresult;
    }
    else
    {
        // Expand master data
        tmp<Field<Type> > expandmaster = expandMasterData(pf);

        tmp<Field<Type> > tresult = patchToPatch().masterToSlave(expandmaster);

        // Truncate to size
        tresult().setSize(size());

        return tresult;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::overlapGgiPolyPatch::interpolate(const tmp<Field<Type> >& tpf) const
{
    if (master())
    {
        // Expand slave data
        tmp<Field<Type> > expand = expandSlaveData(tpf());

        tmp<Field<Type> > tresult = patchToPatch().slaveToMaster(expand);

        // Truncate to size
        tresult().setSize(size());

        return tresult;
    }
    else
    {
        // Expand master data
        tmp<Field<Type> > expandmaster = expandMasterData(tpf());

        tmp<Field<Type> > tresult = patchToPatch().masterToSlave(expandmaster);

        // Truncate to size
        tresult().setSize(size());

        return tresult;
    }
}


// ************************************************************************* //
