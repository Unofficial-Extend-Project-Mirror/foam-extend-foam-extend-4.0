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
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "mixingPlanePolyPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::mixingPlanePolyPatch::interpolate
(
    const Field<Type>& pf
) const
{
    if (master())
    {
        return patchToPatch().slaveToMaster(pf);
    }
    else
    {
        return patchToPatch().masterToSlave(pf);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::mixingPlanePolyPatch::interpolate
(
    const tmp<Field<Type> >& tpf
) const
{
    tmp<Field<Type> > tint = interpolate(tpf());
    tpf.clear();
    return tint;
}

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::mixingPlanePolyPatch::circumferentialAverage
(
    const Field<Type>& pf
) const
{
    if (master())
    {
        return patchToPatch().masterToMaster(pf);
    }
    else
    {
        return patchToPatch().slaveToSlave(pf);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::mixingPlanePolyPatch::circumferentialAverage
(
    const tmp<Field<Type> >& tpf
) const
{
    tmp<Field<Type> > tint = circumferentialAverage(tpf());
    tpf.clear();
    return tint;
}


// ************************************************************************* //
