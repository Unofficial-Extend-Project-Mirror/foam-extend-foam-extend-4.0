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
    Foam::NoHeatTransfer

Description
    Dummy heat transfer model for 'no heat transfer'

\*---------------------------------------------------------------------------*/

#ifndef NoHeatTransfer_H
#define NoHeatTransfer_H

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class NoHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoHeatTransfer
:
    public HeatTransferModel<CloudType>
{

public:

    //- Runtime type information
    TypeName("NoHeatTransfer");


    // Constructors

        //- Construct from dictionary
        NoHeatTransfer
        (
            const dictionary& dict,
            CloudType& cloud
        );


    // Destructor

        ~NoHeatTransfer();


    // Member Functions

        //- Flag to indicate whether model activates heat transfer model
        bool active() const;

        //- Nusselt number
        scalar Nu
        (
            const scalar,
            const scalar
        ) const;

        //- Prandtl number
        scalar Pr() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoHeatTransfer.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
