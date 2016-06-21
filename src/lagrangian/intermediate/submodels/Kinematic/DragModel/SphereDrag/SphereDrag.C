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

\*---------------------------------------------------------------------------*/

#include "SphereDrag.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::SphereDrag<CloudType>::SphereDrag
(
    const dictionary& dict,
    CloudType& owner
)
:
    DragModel<CloudType>(dict, owner)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::SphereDrag<CloudType>::~SphereDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
bool Foam::SphereDrag<CloudType>::active() const
{
    return true;
}


template <class CloudType>
Foam::scalar Foam::SphereDrag<CloudType>::Cd(const scalar Re) const
{
    scalar Cd;
    if (Re < SMALL)
    {
        Cd = GREAT;
    }
    else if (Re > 1000.0)
    {
        Cd =  0.424;
    }
    else
    {
        Cd = 24.0/Re*(1.0 + 1.0/6.0*pow(Re, 2.0/3.0));
    }

    return Cd;
}


// ************************************************************************* //
