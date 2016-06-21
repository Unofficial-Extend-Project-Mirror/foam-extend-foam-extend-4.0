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
    Foam::standardDragModel

Description
    standard drag law model with a lower limit on the Cd-value

\*---------------------------------------------------------------------------*/

#ifndef standardDragModel_H
#define standardDragModel_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class standardDragModel Declaration
\*---------------------------------------------------------------------------*/

class standardDragModel
:
    public dragModel
{

private:

    // Private data

        dictionary dragDict_;
        scalar preReFactor_;
        scalar ReExponent_;
        scalar ReLimiter_;
        scalar CdLimiter_;
        scalar Cdistort_;

public:

    //- Runtime type information
        TypeName("standardDragModel");


    // Constructors

        //- Construct from components
        standardDragModel
        (
            const dictionary& dict
        );


    // Destructor

        ~standardDragModel();


    // Member Functions

        scalar Cd
        (
            const scalar ReynoldsNumber,
            const scalar deviation
        ) const;

        scalar relaxationTime
        (
            const vector& URel,
            const scalar diameter,
            const scalar rho,
            const scalar liquidDensity,
            const scalar nu,
            const scalar deviation
        ) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
