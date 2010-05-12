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

#include "DragModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DragModel<CloudType>::DragModel
(
    const dictionary& dict,
    CloudType& owner
)
:   dict_(dict),
    owner_(owner)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DragModel<CloudType>::~DragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType& Foam::DragModel<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary& Foam::DragModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
Foam::scalar Foam::DragModel<CloudType>::Cu
(
    const vector& Ur,
    const scalar d,
    const scalar rhoc,
    const scalar rhop,
    const scalar mu
) const
{
    const scalar magUr = mag(Ur);

    const scalar Re = rhoc*magUr*d/(mu + SMALL);

    const scalar cd = Cd(Re);

    return 3.0*cd*rhoc*magUr/(4.0*d*rhop);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NewDragModel.C"

// ************************************************************************* //

