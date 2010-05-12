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

\*---------------------------------------------------------------------------*/

#include "basicReactingParcel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicReactingParcel, 0);
    defineParticleTypeNameAndDebug(basicReactingParcel, 0);
    defineParcelTypeNameAndDebug(basicReactingParcel, 0);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicReactingParcel::basicReactingParcel
(
    ReactingCloud<basicReactingParcel>& owner,
    const label typeId,
    const vector& position,
    const label celli,
    const scalar d0,
    const vector& U0,
    const scalar nParticle0,
    const scalarField& YGas0,
    const scalarField& YLiquid0,
    const scalarField& YSolid0,
    const scalarField& YMixture0,
    const constantProperties& constProps
)
:
    ReactingParcel<basicReactingParcel>
    (
        owner,
        typeId,
        position,
        celli,
        d0,
        U0,
        nParticle0,
        YGas0,
        YLiquid0,
        YSolid0,
        YMixture0,
        constProps
    )
{}


Foam::basicReactingParcel::basicReactingParcel
(
    const Cloud<basicReactingParcel>& cloud,
    Istream& is,
    bool readFields
)
:
    ReactingParcel<basicReactingParcel>(cloud, is, readFields)
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

Foam::basicReactingParcel::~basicReactingParcel()
{}


// ************************************************************************* //
