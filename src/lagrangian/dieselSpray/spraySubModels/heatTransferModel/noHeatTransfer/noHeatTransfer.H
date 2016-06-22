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

Class
    Foam::noHeatTransfer

Description
    no heat transfer. return a negative value to prevent
    the liquid being heated up.

\*---------------------------------------------------------------------------*/

#ifndef noHeatTransfer_H
#define noHeatTransfer_H

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class noHeatTransfer
:
    public heatTransferModel
{

public:

    //- Runtime type information
        TypeName("off");


    // Constructors

        //- Construct from dictionary
        noHeatTransfer
        (
            const dictionary& dict
        );


    // Destructor

        ~noHeatTransfer();


    // Member Functions

        bool heatTransfer() const;

        //- Nusselt Number
        scalar Nu
        (
            const scalar ReynoldsNumber,
            const scalar PrandtlNumber
        ) const;

        scalar relaxationTime
        (
            const scalar liquidDensity,
            const scalar diameter,
            const scalar liquidcL,
            const scalar kappa,
            const scalar ReynoldsNumber,
            const scalar PrandtlNumber
        ) const;

        scalar fCorrection(const scalar z) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
