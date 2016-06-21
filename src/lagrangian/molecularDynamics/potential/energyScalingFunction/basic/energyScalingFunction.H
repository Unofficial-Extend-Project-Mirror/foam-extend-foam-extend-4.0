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
    Foam::energyScalingFunction

Description

SourceFiles
    energyScalingFunction.C
    newEnergyScalingFunction.C

\*---------------------------------------------------------------------------*/

#ifndef energyScalingFunction_H
#define energyScalingFunction_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "pairPotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class energyScalingFunction Declaration
\*---------------------------------------------------------------------------*/

class energyScalingFunction
{

protected:

    // Protected data

        word name_;

        dictionary energyScalingFunctionProperties_;

        const pairPotential& pairPot_;


    // Private Member Functions

        //- Disallow copy construct
        energyScalingFunction(const energyScalingFunction&);

        //- Disallow default bitwise assignment
        void operator=(const energyScalingFunction&);


public:

    //- Runtime type information
    TypeName("energyScalingFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            energyScalingFunction,
            dictionary,
            (
                const word& name,
                const dictionary& energyScalingFunctionProperties,
                const pairPotential& pairPot
            ),
            (name, energyScalingFunctionProperties, pairPot)
        );


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<energyScalingFunction> New
        (
            const word& name,
            const dictionary& energyScalingFunctionProperties,
            const pairPotential& pairPot
        );


    // Constructors

        //- Construct from components
        energyScalingFunction
        (
            const word& name,
            const dictionary& energyScalingFunctionProperties,
            const pairPotential& pairPot
        );


    // Destructor

        virtual ~energyScalingFunction()
        {}


    // Member Functions

        virtual void scaleEnergy(scalar& e, const scalar r) const = 0;

        const dictionary& energyScalingFunctionProperties() const
        {
            return energyScalingFunctionProperties_;
        }

        //- Read energyScalingFunction dictionary
        virtual bool read(const dictionary& energyScalingFunctionProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
