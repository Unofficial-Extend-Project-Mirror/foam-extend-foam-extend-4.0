/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
typename Foam::DecoupledBlockCoeff<Type>::scalarType&
Foam::DecoupledBlockCoeff<Type>::toScalar()
{
    if (!scalarCoeffPtr_)
    {
        // Debug check: demotion
        if (linearCoeffPtr_)
        {
            FatalErrorIn
            (
                "DecoupledBlockCoeff<Type>::scalarType& "
                "DecoupledBlockCoeff<Type>::toScalar()"
            )   << "Detected demotion to scalar.  Probably an error"
                << abort(FatalError);
        }

        scalarCoeffPtr_ = new scalarType(pTraits<scalarType>::zero);
    }

    return *scalarCoeffPtr_;
}


template<class Type>
typename Foam::DecoupledBlockCoeff<Type>::linearType&
Foam::DecoupledBlockCoeff<Type>::toLinear()
{
    if (!linearCoeffPtr_)
    {
        linearCoeffPtr_ = new linearType(pTraits<linearType>::zero);

        // If scalar is active, promote to linear
        if (scalarCoeffPtr_)
        {
            *linearCoeffPtr_ = (*scalarCoeffPtr_)*pTraits<linearType>::one;
            deleteDemandDrivenData(scalarCoeffPtr_);
        }
    }

    return *linearCoeffPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::DecoupledBlockCoeff<Type>::DecoupledBlockCoeff()
:
    scalarCoeffPtr_(NULL),
    linearCoeffPtr_(NULL)
{}


template<class Type>
Foam::DecoupledBlockCoeff<Type>::DecoupledBlockCoeff
(
    const DecoupledBlockCoeff<Type>& f
)
:
    scalarCoeffPtr_(NULL),
    linearCoeffPtr_(NULL)
{
    if (f.scalarCoeffPtr_)
    {
        scalarCoeffPtr_ = new scalarType(*(f.scalarCoeffPtr_));
    }
    else if (f.linearCoeffPtr_)
    {
        linearCoeffPtr_ = new linearType(*(f.linearCoeffPtr_));
    }
}


template<class Type>
Foam::DecoupledBlockCoeff<Type>::DecoupledBlockCoeff(Istream& is)
:
    scalarCoeffPtr_(NULL),
    linearCoeffPtr_(NULL)
{
    // Read keyword and pick up allocated field
    word key(is);

    if
    (
        key
     == blockCoeffBase::activeLevelNames_[blockCoeffBase::UNALLOCATED]
    )
    {
    }
    else if
    (
        key
     == blockCoeffBase::activeLevelNames_[blockCoeffBase::SCALAR]
    )
    {
        scalarCoeffPtr_ = new scalarType(readScalar(is));
    }
    else if
    (
        key
     == blockCoeffBase::activeLevelNames_[blockCoeffBase::LINEAR]
    )
    {
        linearCoeffPtr_ = new linearType(is);
    }
    else
    {
        FatalIOErrorIn
        (
            "DecoupledBlockCoeff<Type>::DecoupledBlockCoeff(Istream& is)",
            is
        )   << "invalid keyword while reading: " << key
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::DecoupledBlockCoeff<Type> Foam::DecoupledBlockCoeff<Type>::clone() const
{
    return DecoupledBlockCoeff<Type>(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::DecoupledBlockCoeff<Type>::~DecoupledBlockCoeff()
{
    this->clear();
}


template<class Type>
void Foam::DecoupledBlockCoeff<Type>::clear()
{
    deleteDemandDrivenData(scalarCoeffPtr_);
    deleteDemandDrivenData(linearCoeffPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::blockCoeffBase::activeLevel
Foam::DecoupledBlockCoeff<Type>::activeType() const
{
    if (scalarCoeffPtr_)
    {
        return blockCoeffBase::SCALAR;
    }
    else if (linearCoeffPtr_)
    {
        return blockCoeffBase::LINEAR;
    }
    else
    {
        return blockCoeffBase::UNALLOCATED;
    }
}


template<class Type>
void Foam::DecoupledBlockCoeff<Type>::checkActive() const
{
    label nActive = 0;

    if (scalarCoeffPtr_) nActive++;
    if (linearCoeffPtr_) nActive++;

    if (nActive > 1)
    {
        FatalErrorIn
        (
            "void Foam::DecoupledBlockCoeff<Type>::checkActive() const"
        )   << "Activation/deactivation error.  nActive = " << nActive
            << abort(FatalError);
    }
}


template<class Type>
const typename Foam::DecoupledBlockCoeff<Type>::scalarType&
Foam::DecoupledBlockCoeff<Type>::asScalar() const
{
    if (!scalarCoeffPtr_)
    {
        FatalErrorIn
        (
            "DecoupledBlockCoeff<Type>::scalarType& "
            "DecoupledBlockCoeff<Type>::asScalar()"
        )   << "Requested scalar but active type is: "
            << blockCoeffBase::activeLevelNames_[this->activeType()]
            << ".  This is not allowed."
            << abort(FatalError);
    }

    return *scalarCoeffPtr_;
}


template<class Type>
const typename Foam::DecoupledBlockCoeff<Type>::linearType&
Foam::DecoupledBlockCoeff<Type>::asLinear() const
{
    if (!linearCoeffPtr_)
    {
        FatalErrorIn
        (
            "DecoupledBlockCoeff<Type>::linearType& "
            "DecoupledBlockCoeff<Type>::asLinear()"
        )   << "Requested linear but active type is: "
            << blockCoeffBase::activeLevelNames_[this->activeType()]
            << ".  This is not allowed."
            << abort(FatalError);
    }

    return *linearCoeffPtr_;
}


template<class Type>
typename Foam::DecoupledBlockCoeff<Type>::scalarType&
Foam::DecoupledBlockCoeff<Type>::asScalar()
{
    if (linearCoeffPtr_)
    {
        FatalErrorIn
        (
            "DecoupledBlockCoeff<Type>::scalarType& "
            "DecoupledBlockCoeff<Type>::asScalar()"
        )   << "Requested scalar but active type is: "
            << blockCoeffBase::activeLevelNames_[this->activeType()]
            << ".  This is not allowed."
            << abort(FatalError);
    }

    if (!scalarCoeffPtr_)
    {
        return this->toScalar();
    }

    return *scalarCoeffPtr_;
}


template<class Type>
typename Foam::DecoupledBlockCoeff<Type>::linearType&
Foam::DecoupledBlockCoeff<Type>::asLinear()
{
    if (!linearCoeffPtr_)
    {
        return this->toLinear();
    }

    return *linearCoeffPtr_;
}


template<class Type>
typename Foam::DecoupledBlockCoeff<Type>::scalarType
Foam::DecoupledBlockCoeff<Type>::component(const direction dir) const
{
    if (scalarCoeffPtr_)
    {
        return *scalarCoeffPtr_;
    }
    else if (linearCoeffPtr_)
    {
        return linearCoeffPtr_->component(dir);
    }
    else
    {
        FatalErrorIn
        (
            "tmp<DecoupledBlockCoeff<Type>::scalarType>"
            "DecoupledBlockCoeff<Type>::component(const direction dir) const"
        )   << " not allocated."
            << abort(FatalError);
    }

    // Dummy return to keep compiler happy
    return *scalarCoeffPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::DecoupledBlockCoeff<Type>::operator=
(
    const DecoupledBlockCoeff<Type>& f
)
{
    if (this == &f)
    {
        FatalErrorIn
        (
            "DecoupledBlockCoeff<Type>::operator=("
            "const DecoupledBlockCoeff<Type>&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    if (f.scalarCoeffPtr_)
    {
        this->toScalar() = *(f.scalarCoeffPtr_);
    }
    else if (f.linearCoeffPtr_)
    {
        this->toLinear() = *(f.linearCoeffPtr_);
    }
    else
    {
        // Not allocated - do nothing
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const DecoupledBlockCoeff<Type>& f)
{
    // Write active type
    os << blockCoeffBase::activeLevelNames_[f.activeType()] << nl;

    if (f.scalarCoeffPtr_)
    {
        os << *(f.scalarCoeffPtr_);
    }
    else if (f.linearCoeffPtr_)
    {
        os << *(f.linearCoeffPtr_);
    }

    return os;
}


// ************************************************************************* //
