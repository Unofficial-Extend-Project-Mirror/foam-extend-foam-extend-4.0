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
    Foam::compressibilityModels::linear

Description
    linear compressibility model.

SourceFiles
    linear.C

\*---------------------------------------------------------------------------*/

#ifndef linear_H
#define linear_H

#include "barotropicCompressibilityModel.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressibilityModels
{

/*---------------------------------------------------------------------------*\
                           Class linear Declaration
\*---------------------------------------------------------------------------*/

class linear
:
    public barotropicCompressibilityModel
{
    // Private data

        dimensionedScalar psiv_;
        dimensionedScalar psil_;


public:

    //- Runtime type information
    TypeName("linear");


    // Constructors

        //- construct from components
        linear
        (
            const dictionary& compressibilityProperties,
            const volScalarField& gamma,
            const word& psiName = "psi"
        );


    // Destructor

        ~linear()
        {}


    // Member Functions

        //- Correct the linear compressibility
        void correct();

        //- Read transportProperties dictionary
        bool read(const dictionary& compressibilityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressibilityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
