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

Class
    Foam::NoSurfaceReaction

Description
    Dummy Devolatisation model for 'no devolatisation'

\*---------------------------------------------------------------------------*/

#ifndef NoSurfaceReaction_H
#define NoSurfaceReaction_H

#include "SurfaceReactionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                            Class NoSurfaceReaction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoSurfaceReaction
:
    public SurfaceReactionModel<CloudType>
{

public:

    //- Runtime type information
    TypeName("NoSurfaceReaction");


    // Constructors

        //- Construct from dictionary
        NoSurfaceReaction
        (
            const dictionary& dict,
            CloudType& owner
        );


    // Destructor

        ~NoSurfaceReaction();


    // Member Functions

        //- Flag to indicate whether model activates devolatisation model
        bool active() const;

        //- Update surface reactions
        void calculate
        (
            const scalar dt,
            const label celli,
            const scalar dp,
            const scalar T0,
            const scalar T1,
            const scalar Tc,
            const scalar rhoc,
            const scalar massp,
            scalarField& YGas,
            scalarField& YLiquid,
            scalarField& YSolid,
            scalarField& YMixture,
            scalar& dMassMTSR,
            scalarList& dMassSR
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoSurfaceReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
