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
    Foam::CompositionModel


Description
    Templated reacting parcel composition model class
    Consists of gases (via thermo package), liquids and solids

SourceFiles
    CompositionModel.C
    NewCompositionModel.C

\*---------------------------------------------------------------------------*/

#ifndef CompositionModel_H
#define CompositionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "PtrList.H"
#include "hCombustionThermo.H"

#include "liquidMixture.H"
#include "solidMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class CompositionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CompositionModel
{

    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        //- Reference to the owner injection class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary& coeffDict_;

        //- Reference to the carrier phase thermo package
        hCombustionThermo& carrierThermo_;

        //- Reference to the carrier phase species
        const PtrList<specieReactingProperties>& gases_;

        //- Global liquid properties data
        autoPtr<liquidMixture> liquids_;

        //- Global solid properties data
        autoPtr<solidMixture> solids_;


public:

    //- Runtime type information
    TypeName("CompositionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CompositionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct from dictionary
        CompositionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    // Destructor

        virtual ~CompositionModel();


    // Selector

        static autoPtr<CompositionModel<CloudType> > New
        (
            const dictionary& dict,
            CloudType& owner
        );


    // Member Functions

        // Access

            //- Return the cloud object
            const CloudType& owner() const;

            //- Return the cloud dictionary
            const dictionary& dict() const;

            //- Return the coefficients dictionary
            const dictionary& coeffDict() const;

            //- Return the carrier phase thermo package
            const hCombustionThermo& carrierThermo() const;

            //- Return the gases
            const PtrList<specieReactingProperties>& gases() const;

            //- Return the global liquids
            const liquidMixture& liquids() const;

            //- Return the global solids
            const solidMixture& solids() const;

            //- Return the list of composition names
            virtual const wordList compositionNames() const = 0;

            //- Return the list of gas names
            virtual const wordList& gasNames() const = 0;

            //- Return local id of gas based on gasName
            virtual label gasLocalId(const word& gasName) const = 0;

            //- Return global id of gas based on gasName
            virtual label gasGlobalId(const word& gasName) const = 0;

            //- Return the list indices of gases in global thermo list
            virtual const labelList& gasGlobalIds() const = 0;

            //- Return the list of gas mass fractions
            virtual const scalarField& YGas0() const = 0;

            //- Return the list of liquid names
            virtual const wordList& liquidNames() const = 0;

            //- Return local id of liquid based on liquidName
            virtual label liquidLocalId(const word& liquidName) const = 0;

            //- Return global id of liquid based on liquidName
            virtual label liquidGlobalId(const word& liquidName) const = 0;

            //- Return the list indices of liquid in global thermo list
            virtual const labelList& liquidGlobalIds() const = 0;

            //- Return the list of liquid mass fractions
            virtual const scalarField& YLiquid0() const = 0;

            //- Return the list of solid names
            virtual const wordList& solidNames() const = 0;

            //- Return local id of solid based on solidName
            virtual label solidLocalId(const word& solidName) const = 0;

            //- Return global id of solid based on solidName
            virtual label solidGlobalId(const word& solidName) const = 0;

            //- Return the list indices of solids in global thermo list
            virtual const labelList& solidGlobalIds() const = 0;

            //- Return the list of solid mass fractions
            virtual const scalarField& YSolid0() const = 0;

            //- Return the list of mixture mass fractions
            virtual const scalarField& YMixture0() const = 0;

            //- Return the gas constant for the gas mixture
            virtual scalar RGas(const scalarField& YGas) const = 0;

            //- Return enthalpy for the gas mixture
            virtual scalar HGas
            (
                 const scalarField& YGas,
                 const scalar T
            ) const = 0;

            //- Return specific heat caparcity for the gas mixture
            virtual scalar cpGas
            (
                 const scalarField& YGas,
                 const scalar T
            ) const = 0;

            //- Return specific heat caparcity for the liquid mixture
            virtual scalar cpLiquid
            (
                const scalarField& YLiquid,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Return specific heat caparcity for the solid mixture
            virtual scalar cpSolid
            (
                const scalarField& YSolid
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCompositionModel(CloudType)                                       \
                                                                              \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        CompositionModel<CloudType>,                                          \
        0                                                                     \
    );                                                                        \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        CompositionModel<CloudType>,                                          \
        dictionary                                                            \
    );


#define makeCompositionModelType(SS, CloudType, ParcelType)                   \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType<ParcelType> >, 0);       \
                                                                              \
    CompositionModel<CloudType<ParcelType> >::                                \
        adddictionaryConstructorToTable<SS<CloudType<ParcelType> > >          \
            add##SS##CloudType##ParcelType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CompositionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
