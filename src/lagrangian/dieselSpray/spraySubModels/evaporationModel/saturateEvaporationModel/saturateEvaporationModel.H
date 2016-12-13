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
    Foam::saturateEvaporationModel

Description
    saturate evaporation model.
    Characteristic time is calculated to immediately saturate the cell
    It should be mentioned that this is coupled with the (implicit) way
    evaporation is done in parcel, since the evaporation time depends
    on the integration step

\*---------------------------------------------------------------------------*/

#ifndef saturateEvaporationModel_H
#define saturateEvaporationModel_H

#include "evaporationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class saturateEvaporationModel Declaration
\*---------------------------------------------------------------------------*/

class saturateEvaporationModel
:
    public evaporationModel
{

private:

    // Private data

        dictionary evapDict_;

public:

    //- Runtime type information
        TypeName("saturateEvaporationModel");


    // Constructors

        //- Construct from dictionary
        saturateEvaporationModel
        (
            const dictionary& dict
        );


    // Destructor

        ~saturateEvaporationModel();


    // Member Functions

        bool evaporation() const;

       //- Correlation for the Sherwood Number
        scalar Sh
        (
            const scalar ReynoldsNumber,
            const scalar SchmidtNumber
        ) const;

        //- Return the evaporation relaxation time
        scalar relaxationTime
        (
            const scalar diameter,
            const scalar liquidDensity,
            const scalar rhoFuelVapor,
            const scalar massDiffusionCoefficient,
            const scalar ReynoldsNumber,
            const scalar SchmidtNumber,
            const scalar Xs,
            const scalar Xf,
            const scalar m0,
            const scalar dm,
            const scalar dt
        ) const;

        scalar boilingTime
        (
            const scalar liquidDensity,
            const scalar cpFuel,
            const scalar heatOfVapour,
            const scalar kappa,
            const scalar Nusselt,
            const scalar deltaTemp,
            const scalar diameter,
            const scalar,
            const scalar,
            const scalar,
            const scalar,
            const scalar,
            const scalar,
            const scalar,
            const scalar,
            const scalar
        ) const;

        inline label nEvapIter() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
