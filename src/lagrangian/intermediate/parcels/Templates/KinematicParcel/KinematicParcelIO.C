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

#include "KinematicParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ParcelType>
Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const Cloud<ParcelType>& cloud,
    Istream& is,
    bool readFields
)
:
    Particle<ParcelType>(cloud, is, readFields),
    typeId_(0),
    d_(0.0),
    U_(vector::zero),
    nParticle_(0.0),
    rho_(0.0),
    tTurb_(0.0),
    UTurb_(vector::zero),
    rhoc_(0.0),
    Uc_(vector::zero),
    muc_(0.0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            typeId_ = readLabel(is);
            d_ = readScalar(is);
            is >> U_;
            nParticle_ = readScalar(is);
            rho_ = readScalar(is);
            tTurb_ = readScalar(is);
            is >> UTurb_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&typeId_),
                sizeof(typeId_)
              + sizeof(d_)
              + sizeof(U_)
              + sizeof(nParticle_)
              + sizeof(rho_)
              + sizeof(tTurb_)
              + sizeof(UTurb_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "KinematicParcel<ParcelType>::KinematicParcel"
        "(const Cloud<ParcelType>&, Istream&, bool)"
    );
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::readFields
(
    KinematicCloud<ParcelType>& c
)
{
    if (!c.size())
    {
        return;
    }

    IOField<label> typeId(c.fieldIOobject("typeId"));
    c.checkFieldIOobject(c, typeId);

    IOField<scalar> d(c.fieldIOobject("d"));
    c.checkFieldIOobject(c, d);

    IOField<vector> U(c.fieldIOobject("U"));
    c.checkFieldIOobject(c, U);

    IOField<scalar> nParticle(c.fieldIOobject("nParticle"));
    c.checkFieldIOobject(c, nParticle);

    IOField<scalar> rho(c.fieldIOobject("rho"));
    c.checkFieldIOobject(c, rho);

    IOField<scalar> tTurb(c.fieldIOobject("tTurb"));
    c.checkFieldIOobject(c, tTurb);

    IOField<vector> UTurb(c.fieldIOobject("UTurb"));
    c.checkFieldIOobject(c, UTurb);

    label i = 0;
    forAllIter(typename Cloud<ParcelType>, c, iter)
    {
        ParcelType& p = iter();

        p.typeId_ = typeId[i];
        p.d_ = d[i];
        p.U_ = U[i];
        p.nParticle_ = nParticle[i];
        p.rho_ = rho[i];
        p.tTurb_ = tTurb[i];
        p.UTurb_ = UTurb[i];
        i++;
    }
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::writeFields
(
    const KinematicCloud<ParcelType>& c
)
{
    Particle<ParcelType>::writeFields(c);

    label np =  c.size();

    IOField<label> typeId(c.fieldIOobject("typeId"), np);
    IOField<scalar> d(c.fieldIOobject("d"), np);
    IOField<vector> U(c.fieldIOobject("U"), np);
    IOField<scalar> nParticle(c.fieldIOobject("nParticle"), np);
    IOField<scalar> rho(c.fieldIOobject("rho"), np);
    IOField<scalar> tTurb(c.fieldIOobject("tTurb"), np);
    IOField<vector> UTurb(c.fieldIOobject("UTurb"), np);

    label i = 0;
    forAllConstIter(typename Cloud<ParcelType>, c, iter)
    {
        const KinematicParcel<ParcelType>& p = iter();

        typeId[i] = p.typeId();
        d[i] = p.d();
        U[i] = p.U();
        nParticle[i] = p.nParticle();
        rho[i] = p.rho();
        tTurb[i] = p.tTurb();
        UTurb[i] = p.UTurb();
        i++;
    }

    typeId.write();
    d.write();
    U.write();
    nParticle.write();
    rho.write();
    tTurb.write();
    UTurb.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const KinematicParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const Particle<ParcelType>& >(p)
            << token::SPACE << p.typeId()
            << token::SPACE << p.d()
            << token::SPACE << p.U()
            << token::SPACE << p.nParticle()
            << token::SPACE << p.rho()
            << token::SPACE << p.tTurb()
            << token::SPACE << p.UTurb();
    }
    else
    {
        os  << static_cast<const Particle<ParcelType>& >(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.typeId_),
            sizeof(p.typeId())
          + sizeof(p.d())
          + sizeof(p.U())
          + sizeof(p.nParticle())
          + sizeof(p.rho())
          + sizeof(p.tTurb())
          + sizeof(p.UTurb())
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const KinematicParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
