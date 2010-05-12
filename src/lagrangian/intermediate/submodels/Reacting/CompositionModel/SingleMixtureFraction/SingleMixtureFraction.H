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
    Foam::SingleMixtureFraction


Description
    Templated parcel single mixture fraction class
    - Each phase sums to a mass fraction of 1

SourceFiles
    SingleMixtureFraction.C

\*---------------------------------------------------------------------------*/

#ifndef SingleMixtureFraction_H
#define SingleMixtureFraction_H

#include "CompositionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class SingleMixtureFraction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SingleMixtureFraction
:
    public CompositionModel<CloudType>
{

    // Private data

        //- Parcel properties

            //- List of gas names
            const wordList gasNames_;

            //- Return local id of gas based on gasName
            label gasLocalId(const word& gasName) const;

            //- Return global id of gas based on gasName
            label gasGlobalId(const word& gasName) const;

            //- List of IDs of gases in global list
            labelList gasGlobalIds_;

            //- Initial mass fractions of gases
            const scalarField YGas0_;

            //- Initial total mass fractions of gases
            const scalar YGasTot0_;

            //- List of liquid names
            const wordList liquidNames_;

            //- Return local id of liquid based on liquidName
            label liquidLocalId(const word& liquidName) const;

            //- Return global id of liquid based on liquidName
            label liquidGlobalId(const word& liquidName) const;

            //- List of IDs of liquids in global list
            labelList liquidGlobalIds_;

            //- Initial mass fractions of liquids
            const scalarField YLiquid0_;

            //- Initial total mass fractions of liquids
            const scalar YLiquidTot0_;

            //- List of solids species
            const wordList solidNames_;

            //- Return local id of solid based on solidName
            label solidLocalId(const word& solidName) const;

            //- Return global id of solid based on solidName
            label solidGlobalId(const word& solidName) const;

            //- List of IDs of solids in global list
            labelList solidGlobalIds_;

            //- Initial mass fractions of solids
            const scalarField YSolid0_;

            //- Initial total mass fractions of solids
            const scalar YSolidTot0_;

            //- Initial total mass fractions of the mixture
            //  Index 0 = gas, 1 = liquid, 2 = solid
            scalarField YMixture0_;


public:

    //- Runtime type information
    TypeName("SingleMixtureFraction");


    // Constructors

        //- Construct from dictionary
        SingleMixtureFraction
        (
            const dictionary& dict,
            CloudType& owner
        );


    // Destructor

        virtual ~SingleMixtureFraction();


    // Member Functions

        // Access

            //- Return the list of composition names
            const wordList compositionNames() const;

            //- Return the list of gas names
            const wordList& gasNames() const;

            //- Return the list indices of gases in global thermo list
            const labelList& gasGlobalIds() const;

            //- Return the list of gas mass fractions
            const scalarField& YGas0() const;

            //- Return the total gas mass fraction
            scalar YGasTot0() const;

            //- Return the list of liquid names
            const wordList& liquidNames() const;

            //- Return the list indices of liquid in global thermo list
            const labelList& liquidGlobalIds() const;

            //- Return the list of liquid mass fractions
            const scalarField& YLiquid0() const;

            //- Return the total liquid mass fraction
            scalar YLiquidTot0() const;

             //- Return the list of solid names
            const wordList& solidNames() const;

            //- Return the list indices of solids in global thermo list
            const labelList& solidGlobalIds() const;

            //- Return the list of solid mass fractions
            const scalarField& YSolid0() const;

            //- Return the total solid mass fraction
            scalar YSolidTot0() const;

            //- Return the list of mixture mass fractions
            const scalarField& YMixture0() const;

            //- Return the gas constant for the gas mixture
            scalar RGas(const scalarField& YGas) const;

            //- Return enthalpy for the gas mixture
            scalar HGas(const scalarField& YGas, const scalar T) const;

            //- Return specific heat caparcity for the gas mixture
            scalar cpGas(const scalarField& YGas, const scalar T) const;

            //- Return specific heat caparcity for the liquid mixture
            scalar cpLiquid
            (
                const scalarField& YLiquid,
                const scalar p,
                const scalar T
            ) const;

            //- Return specific heat caparcity for the solid mixture
            scalar cpSolid(const scalarField& YSolid) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SingleMixtureFraction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
