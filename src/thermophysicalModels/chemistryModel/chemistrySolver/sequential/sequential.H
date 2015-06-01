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

Class
    Foam::sequential

Description
    Foam::sequential

SourceFiles
    sequentialI.H
    sequential.C
    sequentialIO.C

\*---------------------------------------------------------------------------*/

#ifndef sequential_H
#define sequential_H

#include "chemistrySolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class CompType, class ThermoType>
class sequential;

/*---------------------------------------------------------------------------*\
                         Class sequential Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class sequential
:
    public chemistrySolver<CompType, ThermoType>
{
    // Private data

        dictionary coeffsDict_;

        // Model constants

            scalar cTauChem_;
            Switch equil_;


public:

    //- Runtime type information
    TypeName("sequential");


    // Constructors


        //- Construct from components
        sequential
        (
            ODEChemistryModel<CompType, ThermoType>& model,
            const word& modelName
        );


    //- Destructor
    virtual ~sequential();


    // Member Functions

        //- Update the concentrations and return the chemical time
        scalar solve
        (
            scalarField &c,
            const scalar T,
            const scalar p,
            const scalar t0,
            const scalar dt
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "sequential.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
