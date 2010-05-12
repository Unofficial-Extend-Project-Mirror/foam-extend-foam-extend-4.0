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

InClass
    Foam::makeReactionJanaf

Description

\*---------------------------------------------------------------------------*/

#ifndef makeReactionJanaf_H
#define makeReactionJanaf_H

#include "Reaction.H"

#include "IrreversibleReaction.H"
#include "ReversibleReaction.H"
#include "NonEquilibriumReversibleReaction.H"

#include "sutherlandTransport.H"
#include "specieThermo.H"
#include "janafThermo.H"
#include "perfectGas.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef sutherlandTransport<specieThermo<janafThermo<perfectGas> > >
    reactionThermo;

typedef Reaction<reactionThermo>
    reaction;

#define makeReaction(Reaction, ReactionRate)                                  \
                                                                              \
    typedef Reaction<reactionThermo, ReactionRate>                            \
        Reaction##ReactionRate;                                               \
                                                                              \
    template<>                                                                \
    const word Reaction##ReactionRate::typeName                               \
    (                                                                         \
        Reaction::typeName_()                                                 \
      + ReactionRate::type()                                                  \
      + reaction::typeName_()                                                 \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable(reaction, Reaction##ReactionRate, Istream);

#define makePressureDependentReaction(Reaction, PressureDependentReactionRate, ReactionRate, FallOffFunction)  \
                                                                              \
    typedef PressureDependentReactionRate<ReactionRate, FallOffFunction>      \
    PressureDependentReactionRate##ReactionRate##FallOffFunction;             \
                                                                              \
    makeReaction                                                              \
    (Reaction, PressureDependentReactionRate##ReactionRate##FallOffFunction)

#define makeIRReactions(ReactionRate)                                         \
    makeReaction(IrreversibleReaction, ReactionRate)                          \
    makeReaction(ReversibleReaction, ReactionRate)

#define makeIRNReactions(ReactionRate)                                        \
    makeIRReactions(ReactionRate)                                             \
    makeReaction(NonEquilibriumReversibleReaction, ReactionRate)

#define makePressureDependentReactions(ReactionRate, FallOffFunction)         \
                                                                              \
    makePressureDependentReaction                                             \
        (IrreversibleReaction, FallOffReactionRate,                           \
         ReactionRate, FallOffFunction)                                       \
                                                                              \
    makePressureDependentReaction                                             \
        (ReversibleReaction, FallOffReactionRate,                             \
         ReactionRate, FallOffFunction)                                       \
                                                                              \
    makePressureDependentReaction                                             \
        (IrreversibleReaction, ChemicallyActivatedReactionRate,               \
         ReactionRate, FallOffFunction)                                       \
                                                                              \
    makePressureDependentReaction                                             \
        (ReversibleReaction, ChemicallyActivatedReactionRate,                 \
         ReactionRate, FallOffFunction)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
