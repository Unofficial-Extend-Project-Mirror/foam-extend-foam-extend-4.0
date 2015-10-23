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

InClass
    Foam::makeReactionThermo

Description
    Macros for instantiating reactions on given thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeReactionThermo_H
#define makeReactionThermo_H

#include "Reaction.H"

#include "IrreversibleReaction.H"
#include "ReversibleReaction.H"
#include "NonEquilibriumReversibleReaction.H"

#include "specieThermo.H"

#include "sutherlandTransport.H"
#include "janafThermo.H"
#include "perfectGas.H"

#include "polynomialTransport.H"
#include "hPolynomialThermo.H"
#include "icoPolynomial.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReaction(Thermo, ReactionType, ReactionRate)                      \
                                                                              \
    typedef Reaction<Thermo> Reaction##Thermo;                                \
                                                                              \
    typedef ReactionType<Thermo, ReactionRate>                                \
        ReactionType##Thermo##ReactionRate;                                   \
                                                                              \
    template<>                                                                \
    const word ReactionType##Thermo##ReactionRate::typeName                   \
    (                                                                         \
        ReactionType::typeName_()                                             \
      + ReactionRate::type()                                                  \
      + Reaction##Thermo::typeName_()                                 \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        Reaction##Thermo,                                                     \
        ReactionType##Thermo##ReactionRate,                                   \
        Istream                                                               \
    );


#define makePressureDependentReaction(Thermo, Reaction, PressureDependentReactionRate, ReactionRate, FallOffFunction) \
                                                                              \
    typedef PressureDependentReactionRate<ReactionRate, FallOffFunction>      \
        PressureDependentReactionRate##ReactionRate##FallOffFunction;         \
                                                                              \
    makeReaction                                                              \
    (                                                                         \
        Thermo,                                                               \
        Reaction,                                                             \
        PressureDependentReactionRate##ReactionRate##FallOffFunction          \
    )


#define makeIRReactions(Thermo, ReactionRate)                                 \
                                                                              \
    makeReaction(Thermo, IrreversibleReaction, ReactionRate)                  \
                                                                              \
    makeReaction(Thermo, ReversibleReaction, ReactionRate)


#define makeIRNReactions(Thermo, ReactionRate)                                \
                                                                              \
    makeIRReactions(Thermo, ReactionRate)                                     \
                                                                              \
    makeReaction(Thermo, NonEquilibriumReversibleReaction, ReactionRate)


#define makePressureDependentReactions(Thermo, ReactionRate, FallOffFunction) \
                                                                              \
    makePressureDependentReaction                                             \
    (                                                                         \
        Thermo,                                                               \
        IrreversibleReaction,                                                 \
        FallOffReactionRate,                                                  \
        ReactionRate,                                                         \
        FallOffFunction                                                       \
    )                                                                         \
                                                                              \
    makePressureDependentReaction                                             \
    (                                                                         \
        Thermo,                                                               \
        ReversibleReaction,                                                   \
        FallOffReactionRate,                                                  \
        ReactionRate,                                                         \
        FallOffFunction                                                       \
    )                                                                         \
                                                                              \
    makePressureDependentReaction                                             \
    (                                                                         \
        Thermo,                                                               \
        IrreversibleReaction,                                                 \
        ChemicallyActivatedReactionRate,                                      \
        ReactionRate,                                                         \
        FallOffFunction                                                       \
    )                                                                         \
                                                                              \
    makePressureDependentReaction                                             \
    (                                                                         \
        Thermo,                                                               \
        ReversibleReaction,                                                   \
        ChemicallyActivatedReactionRate,                                      \
        ReactionRate,                                                         \
        FallOffFunction                                                       \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
