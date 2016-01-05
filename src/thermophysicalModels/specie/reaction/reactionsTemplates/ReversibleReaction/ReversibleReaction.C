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

#include "ReversibleReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class ReactionThermo, class ReactionRate>
ReversibleReaction<ReactionThermo, ReactionRate>::ReversibleReaction
(
    const Reaction<ReactionThermo>& reaction,
    const ReactionRate& k
)
:
    Reaction<ReactionThermo>(reaction),
    k_(k)
{}


// Construct from components
template<class ReactionThermo, class ReactionRate>
ReversibleReaction<ReactionThermo, ReactionRate>::ReversibleReaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    Istream& is
)
:
    Reaction<ReactionThermo>(species, thermoDatabase, is),
    k_(species, is)
{}


// Construct as copy given new speciesTable
template<class ReactionThermo, class ReactionRate>
ReversibleReaction<ReactionThermo, ReactionRate>::ReversibleReaction
(
    const ReversibleReaction<ReactionThermo, ReactionRate>& rr,
    const speciesTable& species
)
:
    Reaction<ReactionThermo>(rr, species),
    k_(rr.k_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ReactionRate>
scalar ReversibleReaction<ReactionThermo, ReactionRate>::kf
(
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return k_(T, p, c);
}


template<class ReactionThermo, class ReactionRate>
scalar ReversibleReaction<ReactionThermo, ReactionRate>::kr
(
    const scalar kfwd,
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return kfwd/this->Kc(T);
}


template<class ReactionThermo, class ReactionRate>
scalar ReversibleReaction<ReactionThermo, ReactionRate>::kr
(
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return kr(kf(T, p, c), T, p, c);
}


template<class ReactionThermo, class ReactionRate>
void ReversibleReaction<ReactionThermo, ReactionRate>::write
(
    Ostream& os
) const
{
    Reaction<ReactionThermo>::write(os);
    os  << token::SPACE << k_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
