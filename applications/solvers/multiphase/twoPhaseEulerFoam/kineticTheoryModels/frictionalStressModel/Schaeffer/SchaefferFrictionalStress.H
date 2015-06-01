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
    Foam::SchaefferFrictionalStress

Description

SourceFiles
    SchaefferFrictionalStress.C

\*---------------------------------------------------------------------------*/

#ifndef SchaefferFrictionalStress_H
#define SchaefferFrictionalStress_H

#include "frictionalStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SchaefferFrictionalStress Declaration
\*---------------------------------------------------------------------------*/

class SchaefferFrictionalStress
:
    public frictionalStressModel
{

public:

        //- Runtime type information
        TypeName("Schaeffer");


    // Constructors

        //- Construct from components
        SchaefferFrictionalStress(const dictionary& dict);


    //- Destructor
    virtual ~SchaefferFrictionalStress();


    // Member functions

        virtual tmp<volScalarField> frictionalPressure
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax,
            const dimensionedScalar& Fr,
            const dimensionedScalar& eta,
            const dimensionedScalar& p
        ) const;

        virtual tmp<volScalarField> frictionalPressurePrime
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax,
            const dimensionedScalar& Fr,
            const dimensionedScalar& n,
            const dimensionedScalar& p
        ) const;

        virtual tmp<volScalarField> muf
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax,
            const volScalarField& pf,
            const volTensorField& D,
            const dimensionedScalar& phi
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
