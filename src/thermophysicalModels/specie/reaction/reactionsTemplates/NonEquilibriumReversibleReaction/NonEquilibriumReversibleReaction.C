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

#include "NonEquilibriumReversibleReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class ReactionThermo, class ReactionRate>
NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>::
NonEquilibriumReversibleReaction
(
    const Reaction<ReactionThermo>& reaction,
    const ReactionRate& forwardReactionRate,
    const ReactionRate& reverseReactionRate
)
:
    Reaction<ReactionThermo>(reaction),
    fk_(forwardReactionRate),
    rk_(reverseReactionRate)
{}


// Construct from components
template<class ReactionThermo, class ReactionRate>
NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>::
NonEquilibriumReversibleReaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    Istream& is
)
:
    Reaction<ReactionThermo>(species, thermoDatabase, is),
    fk_(species, is),
    rk_(species, is)
{}

// Construct as copy given new speciesTable
template<class ReactionThermo, class ReactionRate>
NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>::
NonEquilibriumReversibleReaction
(
    const NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>& nerr,
    const speciesTable& species
)
:
    Reaction<ReactionThermo>(nerr, species),
    fk_(nerr.fk_),
    rk_(nerr.rk_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ReactionRate>
scalar NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>::kf
(
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return fk_(T, p, c);
}


template<class ReactionThermo, class ReactionRate>
scalar NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>::kr
(
    const scalar,
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return rk_(T, p, c);
}


template<class ReactionThermo, class ReactionRate>
scalar NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>::kr
(
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return rk_(T, p, c);
}


template<class ReactionThermo, class ReactionRate>
void NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>::write
(
    Ostream& os
) const
{
    Reaction<ReactionThermo>::write(os);
    os  << token::SPACE << fk_ << token::SPACE << rk_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
