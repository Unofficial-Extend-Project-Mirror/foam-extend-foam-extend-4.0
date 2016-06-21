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

#include "PhaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
const Foam::wordList Foam::PhaseChangeModel<CloudType>::
enthalpyTransferTypeNames
(
    IStringStream
    (
        "("
            "latentHeat "
            "enthalpyDifference"
        ")"
    )()
);


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class CloudType>
typename Foam::PhaseChangeModel<CloudType>::enthalpyTransferType
Foam::PhaseChangeModel<CloudType>::wordToEnthalpyTransfer(const word& etName)
const
{
    forAll(enthalpyTransferTypeNames, i)
    {
        if (etName == enthalpyTransferTypeNames[i])
        {
            return enthalpyTransferType(i);
        }
    }

    FatalErrorIn
    (
        "PhaseChangeModel<CloudType>::enthalpyTransferType"
        "PhaseChangeModel<CloudType>::wordToEnthalpyTransfer(const word&) const"
    )   << "Unknown enthalpyType " << etName << ". Valid selections are:" << nl
        << enthalpyTransferTypeNames << exit(FatalError);

    return enthalpyTransferType(0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PhaseChangeModel<CloudType>::PhaseChangeModel
(
    CloudType& owner
)
:
    dict_(dictionary::null),
    owner_(owner),
    coeffDict_(dictionary::null),
    enthalpyTransfer_(etLatentHeat)
{}


template<class CloudType>
Foam::PhaseChangeModel<CloudType>::PhaseChangeModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    dict_(dict),
    owner_(owner),
    coeffDict_(dict.subDict(type + "Coeffs")),
    enthalpyTransfer_
    (
        wordToEnthalpyTransfer(coeffDict_.lookup("enthalpyTransfer"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PhaseChangeModel<CloudType>::~PhaseChangeModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class CloudType>
const CloudType& Foam::PhaseChangeModel<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary& Foam::PhaseChangeModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary& Foam::PhaseChangeModel<CloudType>::coeffDict() const
{
    return coeffDict_;
}


template<class CloudType>
const typename Foam::PhaseChangeModel<CloudType>::enthalpyTransferType&
Foam::PhaseChangeModel<CloudType>::enthalpyTransfer() const
{
    return enthalpyTransfer_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NewPhaseChangeModel.C"

// ************************************************************************* //

