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

#include "DispersionRASModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DispersionRASModel<CloudType>::DispersionRASModel
(
    const dictionary& dict,
    CloudType& owner
)
:
    DispersionModel<CloudType>(dict, owner),
    turbulence_
    (
        owner.mesh().objectRegistry::template
        lookupObject<compressible::RASModel>
        (
            "RASProperties"
        )
    ),
    kPtr_(NULL),
    ownK_(false),
    epsilonPtr_(NULL),
    ownEpsilon_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DispersionRASModel<CloudType>::~DispersionRASModel()
{
    cacheFields(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::DispersionRASModel<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        tmp<volScalarField> tk = this->turbulence().k();
        if (tk.isTmp())
        {
            kPtr_ = tk.ptr();
            ownK_ = true;
        }
        else
        {
            kPtr_ = tk.operator->();
            ownK_ = false;
        }

        tmp<volScalarField> tepsilon = this->turbulence().epsilon();
        if (tepsilon.isTmp())
        {
            epsilonPtr_ = tepsilon.ptr();
            ownEpsilon_ = true;
        }
        else
        {
            epsilonPtr_ = tepsilon.operator->();
            ownEpsilon_ = false;
        }
    }
    else
    {
        if (ownK_ && kPtr_)
        {
            delete kPtr_;
            kPtr_ = NULL;
            ownK_ = false;
        }
        if (ownEpsilon_ && epsilonPtr_)
        {
            delete epsilonPtr_;
            epsilonPtr_ = NULL;
            ownEpsilon_ = false;
        }
    }
}


// ************************************************************************* //
