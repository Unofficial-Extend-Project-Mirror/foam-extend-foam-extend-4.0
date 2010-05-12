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

\*---------------------------------------------------------------------------*/

#include "SingleMixtureFraction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SingleMixtureFraction<CloudType>::SingleMixtureFraction
(
    const dictionary& dict,
    CloudType& owner
)
:
    CompositionModel<CloudType>(dict, owner, typeName),

    gasNames_(this->coeffDict().lookup("gasNames")),
    gasGlobalIds_(gasNames_.size(), -1),
    YGas0_(this->coeffDict().lookup("YGas0")),
    YGasTot0_(readScalar(this->coeffDict().lookup("YGasTot0"))),

    liquidNames_(this->coeffDict().lookup("liquidNames")),
    liquidGlobalIds_(liquidNames_.size(), -1),
    YLiquid0_(this->coeffDict().lookup("YLiquid0")),
    YLiquidTot0_(readScalar(this->coeffDict().lookup("YLiquidTot0"))),

    solidNames_(this->coeffDict().lookup("solidNames")),
    solidGlobalIds_(solidNames_.size(), -1),
    YSolid0_(this->coeffDict().lookup("YSolid0")),
    YSolidTot0_(readScalar(this->coeffDict().lookup("YSolidTot0"))),

    YMixture0_(3)
{
    // Construct gasGlobalIds_ list
    forAll(gasNames_, i)
    {
        forAll (this->carrierThermo().composition().Y(), j)
        {
            word specieName(this->carrierThermo().composition().Y()[j].name());

            if (specieName == gasNames_[i])
            {
                gasGlobalIds_[i] = j;
                break;
            }
        }
        if (gasGlobalIds_[i] == -1)
        {
            Info<< "\nThermo package species composition comprises:" << endl;
            forAll (this->carrierThermo().composition().Y(), k)
            {
                Info<< this->carrierThermo().composition().Y()[k].name()
                    << endl;
            }

            FatalErrorIn
            (
                "Foam::SingleMixtureFraction<CloudType>::SingleMixtureFraction"
                "(const dictionary&, CloudType&)"
            )   << "Could not find gas species " << gasNames_[i]
                << " in species list" <<  exit(FatalError);
        }
    }

    // Construct liquidGlobalIds_ list
    forAll(liquidNames_, i)
    {
        forAll (this->liquids().components(), j)
        {
            word specieName(this->liquids().components()[j]);

            if (specieName == liquidNames_[i])
            {
                liquidGlobalIds_[i] = j;
                break;
            }
        }
        if (liquidGlobalIds_[i] == -1)
        {
            Info<< "\nLiquid mixture species composition comprises:" << endl;
            forAll (this->liquids().components(), k)
            {
                Info<< this->liquids().components()[k] << endl;
            }

            FatalErrorIn
            (
                "Foam::SingleMixtureFraction<CloudType>::SingleMixtureFraction"
                "(const dictionary&, CloudType&)"
            )   << "Could not find liquid species " << liquidNames_[i]
                << " in species list" <<  exit(FatalError);
        }
    }

    // Construct solidGlobalIds_ list
    forAll(solidNames_, i)
    {
        forAll (this->solids().components(), j)
        {
            word specieName(this->solids().components()[j]);

            if (specieName == solidNames_[i])
            {
                solidGlobalIds_[i] = j;
                break;
            }
        }
        if (solidGlobalIds_[i] == -1)
        {
            Info<< "\nSolid mixture species composition comprises:" << endl;
            forAll (this->solids().components(), k)
            {
                Info<< this->solids().components()[k] << endl;
            }

            FatalErrorIn
            (
                "Foam::SingleMixtureFraction<CloudType>::SingleMixtureFraction"
                "(const dictionary&, CloudType&)"
            )   << "Could not find solid species " << solidNames_[i]
                << " in species list" <<  exit(FatalError);
        }
    }

    // Set mixture fractions
    YMixture0_[0] = YGasTot0_;
    YMixture0_[1] = YLiquidTot0_;
    YMixture0_[2] = YSolidTot0_;

    // Check that total mass fractions = 1

    if (YGas0_.size())
    {
        if (mag(sum(YGas0_) - 1) > SMALL)
        {
            FatalErrorIn
            (
                "Foam::SingleMixtureFraction<CloudType>::SingleMixtureFraction"
                "(const dictionary&, CloudType&)"
            )<< "Mass fractions of YGas0 must sum to unity"
             <<  exit(FatalError);
        }
    }
    else
    {
        YMixture0_[0] = 0.0;
    }

    if (YLiquid0_.size())
    {
        if (mag(sum(YLiquid0_) - 1) > SMALL)
        {
            FatalErrorIn
            (
                "Foam::SingleMixtureFraction<CloudType>::SingleMixtureFraction"
                "(const dictionary&, CloudType&)"
            )<< "Mass fractions of YLiquid0 must sum to unity"
             <<  exit(FatalError);
        }
    }
    else
    {
        YMixture0_[1] = 0.0;
    }

    if (YSolid0_.size())
    {
        if (mag(sum(YSolid0_) - 1) > SMALL)
        {
            FatalErrorIn
            (
                "Foam::SingleMixtureFraction<CloudType>::SingleMixtureFraction"
                "(const dictionary&, CloudType&)"
            )<< "Mass fractions of YSolid0 must sum to unity"
             <<  exit(FatalError);
        }
    }
    else
    {
        YMixture0_[2] = 0.0;
    }

    // Check total mixture fraction sums to 1
    if (mag(sum(YMixture0_) - 1) > SMALL)
    {
        FatalErrorIn
        (
            "Foam::SingleMixtureFraction<CloudType>::SingleMixtureFraction"
            "(const dictionary&, CloudType&)"
        )   << "Mass fractions YGasTot0 + YSolidTot0 + YSolidTot0 must sum "
            << "to unity" <<  exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SingleMixtureFraction<CloudType>::~SingleMixtureFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const Foam::wordList
Foam::SingleMixtureFraction<CloudType>::compositionNames() const
{
    wordList names(3);
    names[0] = "Gas";
    names[1] = "Liquid";
    names[2] = "Solid";
    return names;
}

template<class CloudType>
const Foam::wordList&
Foam::SingleMixtureFraction<CloudType>::gasNames() const
{
     return gasNames_;
}


template<class CloudType>
Foam::label
Foam::SingleMixtureFraction<CloudType>::gasLocalId(const word& gasName) const
{
    forAll(gasNames_, i)
    {
        if (gasName == gasNames_[i])
        {
            return i;
        }
    }

    WarningIn
    (
        "Foam::label SingleMixtureFraction<CloudType>::"
        "gasLocalId(const word& gasName) const"
    )<< "Gas name " << gasName << " not found in gasNames_"
     << endl;

    return -1;
}


template<class CloudType>
Foam::label
Foam::SingleMixtureFraction<CloudType>::gasGlobalId(const word& gasName) const
{
    forAll(gasNames_, i)
    {
        if (gasName == gasNames_[i])
        {
            return gasGlobalIds_[i];
        }
    }

    WarningIn
    (
        "Foam::label SingleMixtureFraction<CloudType>::"
        "gasGlobalId(const word& gasName) const"
    )<< "Gas name " << gasName << " not found in gasNames_"
     << endl;

    return -1;
}


template<class CloudType>
const Foam::labelList&
Foam::SingleMixtureFraction<CloudType>::gasGlobalIds() const
{
    return gasGlobalIds_;
}


template<class CloudType>
const Foam::scalarField&
Foam::SingleMixtureFraction<CloudType>::YGas0() const
{
    return YGas0_;
}


template<class CloudType>
Foam::scalar Foam::SingleMixtureFraction<CloudType>::YGasTot0() const
{
    return YGasTot0_;
}


template<class CloudType>
const Foam::wordList&
Foam::SingleMixtureFraction<CloudType>::liquidNames() const
{
     return liquidNames_;
}


template<class CloudType>
Foam::label Foam::SingleMixtureFraction<CloudType>::liquidLocalId
(
    const word& liquidName
) const
{
    forAll(liquidNames_, i)
    {
        if (liquidName == liquidNames_[i])
        {
            return i;
        }
    }

    WarningIn
    (
        "Foam::label SingleMixtureFraction<CloudType>::"
        "liquidLocalId(const word& liquidName) const"
    )<< "Liquid name " << liquidName << " not found in liquidNames_"
     << endl;

    return -1;
}


template<class CloudType>
Foam::label Foam::SingleMixtureFraction<CloudType>::liquidGlobalId
(
    const word& liquidName
) const
{
    forAll(liquidNames_, i)
    {
        if (liquidName == liquidNames_[i])
        {
            return liquidGlobalIds_[i];
        }
    }

    WarningIn
    (
        "Foam::label SingleMixtureFraction<CloudType>::"
        "liquidGlobalId(const word& liquidName) const"
    )<< "Liquid name " << liquidName << " not found in liquidNames_"
     << endl;

    return -1;
}


template<class CloudType>
const Foam::labelList&
Foam::SingleMixtureFraction<CloudType>::liquidGlobalIds() const
{
    return liquidGlobalIds_;
}


template<class CloudType>
const Foam::scalarField&
Foam::SingleMixtureFraction<CloudType>::YLiquid0() const
{
    return YLiquid0_;
}


template<class CloudType>
Foam::scalar Foam::SingleMixtureFraction<CloudType>::YLiquidTot0() const
{
    return YLiquidTot0_;
}


template<class CloudType>
const Foam::wordList&
Foam::SingleMixtureFraction<CloudType>::solidNames() const
{
    return solidNames_;
}


template<class CloudType>
Foam::label Foam::SingleMixtureFraction<CloudType>::solidLocalId
(
    const word& solidName
) const
{
    forAll(solidNames_, i)
    {
        if (solidName == solidNames_[i])
        {
            return i;
        }
    }

    WarningIn
    (
        "Foam::label SingleMixtureFraction<CloudType>::"
        "SolididLocalId(const word& solidName) const"
    )<< "Solid name " << solidName << " not found in solidNames_"
     << endl;

    return -1;
}


template<class CloudType>
Foam::label
Foam::SingleMixtureFraction<CloudType>::solidGlobalId
(
    const word& solidName
) const
{
    forAll(solidNames_, i)
    {
        if (solidName == solidNames_[i])
        {
            return solidGlobalIds_[i];
        }
    }

    WarningIn
    (
        "Foam::label SingleMixtureFraction<CloudType>::"
        "solidGlobalId(const word& solidName) const"
    )<< "Solid name " << solidName << " not found in solidNames_"
     << endl;

    return -1;
}


template<class CloudType>
const Foam::labelList&
Foam::SingleMixtureFraction<CloudType>::solidGlobalIds() const
{
    return solidGlobalIds_;
}


template<class CloudType>
const Foam::scalarField&
Foam::SingleMixtureFraction<CloudType>::YSolid0() const
{
    return YSolid0_;
}


template<class CloudType>
Foam::scalar
Foam::SingleMixtureFraction<CloudType>::YSolidTot0() const
{
    return YSolidTot0_;
}


template<class CloudType>
const Foam::scalarField&
Foam::SingleMixtureFraction<CloudType>::YMixture0() const
{
    return YMixture0_;
}


template<class CloudType>
Foam::scalar Foam::SingleMixtureFraction<CloudType>::RGas
(
    const scalarField& YGas
) const
{
    scalar RGasMixture = 0.0;
    forAll(YGas, i)
    {
        label id = gasGlobalIds_[i];
        RGasMixture += YGas[i]*this->gases()[id].R();
    }
    return RGasMixture;
}


template<class CloudType>
Foam::scalar Foam::SingleMixtureFraction<CloudType>::HGas
(
    const scalarField& YGas,
    const scalar T
)
const
{
    scalar HMixture = 0.0;
    forAll(YGas, i)
    {
        label id = gasGlobalIds_[i];
        HMixture += YGas[i]*this->gases()[id].H(T);
    }
    return HMixture;
}


template<class CloudType>
Foam::scalar Foam::SingleMixtureFraction<CloudType>::cpGas
(
    const scalarField& YGas,
    const scalar T
)
const
{
    scalar cpMixture = 0.0;
    forAll(YGas, i)
    {
        label id = gasGlobalIds_[i];
        cpMixture += YGas[i]*this->gases()[id].Cp(T);
    }
    return cpMixture;
}


template<class CloudType>
Foam::scalar Foam::SingleMixtureFraction<CloudType>::cpLiquid
(
    const scalarField& YLiquid,
    const scalar p,
    const scalar T
)
const
{
    scalar cpMixture = 0.0;
    forAll(YLiquid, i)
    {
        label id = liquidGlobalIds_[i];
        cpMixture += YLiquid[i]*this->liquids().properties()[id].cp(p, T);
    }

    return cpMixture;
}


template<class CloudType>
Foam::scalar Foam::SingleMixtureFraction<CloudType>::cpSolid
(
    const scalarField& YSolid
)
const
{
    scalar cpMixture = 0.0;
    forAll(YSolid, i)
    {
        label id = solidGlobalIds_[i];
        cpMixture += YSolid[i]*this->solids().properties()[id].cp();
    }

    return cpMixture;
}


// ************************************************************************* //

