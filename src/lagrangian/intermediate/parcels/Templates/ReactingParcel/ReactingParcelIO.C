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

#include "ReactingParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingParcel<ParcelType>::ReactingParcel
(
    const Cloud<ParcelType>& cloud,
    Istream& is,
    bool readFields
)
:
    ThermoParcel<ParcelType>(cloud, is, readFields),
    mass0_(0.0),
    YMixture_(0),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0),
    pc_(0.0)
{
    if (readFields)
    {
        const ReactingCloud<ParcelType>& cR =
            dynamic_cast<const ReactingCloud<ParcelType>& >(cloud);

        const label nMixture = cR.composition().compositionNames().size();
        const label nGas = cR.composition().gasNames().size();
        const label nLiquid = cR.composition().liquidNames().size();
        const label nSolid = cR.composition().solidNames().size();

        YMixture_.setSize(nMixture);
        YGas_.setSize(nGas);
        YLiquid_.setSize(nLiquid);
        YSolid_.setSize(nSolid);

        if (is.format() == IOstream::ASCII)
        {
            is >> mass0_ >> YMixture_ >> YGas_ >> YLiquid_ >> YSolid_;
            YGas_ /= YMixture_[0] + VSMALL;
            YLiquid_ /= YMixture_[1] + VSMALL;
            YSolid_ /= YMixture_[2] + VSMALL;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&mass0_),
              + sizeof(mass0_)
            );
            is >> YMixture_ >> YGas_ >> YLiquid_ >> YSolid_;
            YGas_ /= YMixture_[0] + VSMALL;
            YLiquid_ /= YMixture_[1] + VSMALL;
            YSolid_ /= YMixture_[2] + VSMALL;
        }
    }

    // Check state of Istream
    is.check
    (
        "ReactingParcel<ParcelType>::ReactingParcel"
        "(const Cloud<ParcelType>&, Istream&, bool)"
    );
}


template<class ParcelType>
void Foam::ReactingParcel<ParcelType>::readFields
(
    ReactingCloud<ParcelType>& c
)
{
    if (!c.size())
    {
        return;
    }

    ThermoParcel<ParcelType>::readFields(c);

    IOField<scalar> mass0(c.fieldIOobject("mass0"));
    c.checkFieldIOobject(c, mass0);

    label i = 0;
    forAllIter(typename Cloud<ParcelType>, c, iter)
    {
        ReactingParcel<ParcelType>& p = iter();
        p.mass0_ = mass0[i++];
    }

    // Get names and sizes for each Y...
    const wordList compositionNames = c.composition().compositionNames();
    const wordList gasNames = c.composition().gasNames();
    const wordList liquidNames = c.composition().liquidNames();
    const wordList solidNames = c.composition().solidNames();
    const label nComposition = compositionNames.size();
    const label nGas = gasNames.size();
    const label nLiquid = liquidNames.size();
    const label nSolid = solidNames.size();

    // Set storage for each Y... for each parcel
    forAllIter(typename Cloud<ParcelType>, c, iter)
    {
        ReactingParcel<ParcelType>& p = iter();
        p.YMixture_.setSize(nComposition, 0.0);
        p.YGas_.setSize(nGas, 0.0);
        p.YLiquid_.setSize(nLiquid, 0.0);
        p.YSolid_.setSize(nSolid, 0.0);
    }

    // Populate YMixture for each parcel
    forAll(compositionNames, j)
    {
        IOField<scalar> YMixture(c.fieldIOobject("Y" + compositionNames[j]));

        label i = 0;
        forAllIter(typename Cloud<ParcelType>, c, iter)
        {
            ReactingParcel<ParcelType>& p = iter();
            p.YMixture_[j] = YMixture[i++];
        }
    }
    // Populate YGas for each parcel
    forAll(gasNames, j)
    {
        IOField<scalar> YGas(c.fieldIOobject("Y" + gasNames[j]));

        label i = 0;
        forAllIter(typename Cloud<ParcelType>, c, iter)
        {
            ReactingParcel<ParcelType>& p = iter();
            p.YGas_[j] = YGas[i++]/p.YMixture_[0];
        }
    }
    // Populate YLiquid for each parcel
    forAll(liquidNames, j)
    {
        IOField<scalar> YLiquid(c.fieldIOobject("Y" + liquidNames[j]));

        label i = 0;
        forAllIter(typename Cloud<ParcelType>, c, iter)
        {
            ReactingParcel<ParcelType>& p = iter();
            p.YLiquid_[j] = YLiquid[i++]/p.YMixture_[1];
        }
    }
    // Populate YSolid for each parcel
    forAll(solidNames, j)
    {
        IOField<scalar> YSolid(c.fieldIOobject("Y" + solidNames[j]));

        label i = 0;
        forAllIter(typename Cloud<ParcelType>, c, iter)
        {
            ReactingParcel<ParcelType>& p = iter();
            p.YSolid_[j] = YSolid[i++]/p.YMixture_[2];
        }
    }
}


template<class ParcelType>
void Foam::ReactingParcel<ParcelType>::writeFields
(
    const ReactingCloud<ParcelType>& c
)
{
    ThermoParcel<ParcelType>::writeFields(c);

    label np =  c.size();

    IOField<scalar> mass0(c.fieldIOobject("mass0"), np);

    label i = 0;
    forAllConstIter(typename Cloud<ParcelType>, c, iter)
    {
        const ReactingParcel<ParcelType>& p = iter();
        mass0[i++] = p.mass0_;
    }

    mass0.write();

    // Write the composition fractions
    if (np > 0)
    {
        const wordList compositionNames = c.composition().compositionNames();
        forAll(compositionNames, j)
        {
            IOField<scalar> YMixture
            (
                c.fieldIOobject("Y" + compositionNames[j]), np
            );

            label i = 0;
            forAllConstIter(typename Cloud<ParcelType>, c, iter)
            {
                const ReactingParcel<ParcelType>& p0 = iter();
                YMixture[i++] = p0.YMixture()[j];
            }

            YMixture.write();
        }
        const wordList& gasNames = c.composition().gasNames();
        forAll(gasNames, j)
        {
            IOField<scalar> YGas(c.fieldIOobject("Y" + gasNames[j]), np);

            label i = 0;
            forAllConstIter(typename Cloud<ParcelType>, c, iter)
            {
                const ReactingParcel<ParcelType>& p0 = iter();
                YGas[i++] = p0.YGas()[j]*p0.YMixture()[0];
            }

            YGas.write();
        }
        const wordList& liquidNames = c.composition().liquidNames();
        forAll(liquidNames, j)
        {
            IOField<scalar> YLiquid(c.fieldIOobject("Y" + liquidNames[j]), np);

            label i = 0;
            forAllConstIter(typename Cloud<ParcelType>, c, iter)
            {
                const ReactingParcel<ParcelType>& p0 = iter();
                YLiquid[i++] = p0.YLiquid()[j]*p0.YMixture()[1];
            }

            YLiquid.write();
        }
        const wordList& solidNames = c.composition().solidNames();
        forAll(solidNames, j)
        {
            IOField<scalar> YSolid(c.fieldIOobject("Y" + solidNames[j]), np);

            label i = 0;
            forAllConstIter(typename Cloud<ParcelType>, c, iter)
            {
                const ReactingParcel<ParcelType>& p0 = iter();
                YSolid[i++] = p0.YSolid()[j]*p0.YMixture()[2];
            }

            YSolid.write();
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ReactingParcel<ParcelType>& p
)
{
    scalarField YGasLoc = p.YGas()*p.YMixture()[0];
    scalarField YLiquidLoc = p.YLiquid()*p.YMixture()[1];
    scalarField YSolidLoc = p.YSolid()*p.YMixture()[2];
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ThermoParcel<ParcelType>& >(p)
            << token::SPACE << p.mass0()
            << token::SPACE << p.YMixture()
            << token::SPACE << YGasLoc
            << token::SPACE << YLiquidLoc
            << token::SPACE << YSolidLoc;
    }
    else
    {
        os  << static_cast<const ThermoParcel<ParcelType>& >(p);
        os.write
        (
            reinterpret_cast<const char*>
            (
                &const_cast<ReactingParcel<ParcelType>&>(p).mass0()
            ),
            sizeof(p.mass0())
        );
        os << p.YMixture() << YGasLoc << YLiquidLoc << YSolidLoc;
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const ReactingParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
