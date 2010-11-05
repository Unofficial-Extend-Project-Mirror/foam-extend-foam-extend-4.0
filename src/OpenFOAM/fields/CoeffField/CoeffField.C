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
#include "expandTensorField.H"

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

template<class Type>
const char* const Foam::CoeffField<Type>::typeName("CoeffField");


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class Type2>
inline void Foam::CoeffField<Type>::checkSize(const UList<Type2>& f) const
{
    if (f.size() != this->size())
    {
        FatalErrorIn
        (
            "void CoeffField<Type>::checkSize(const Field<Type2>& f) const"
        )   << "Incorrect field size: " << f.size()
            << " local size: " << size()
            << abort(FatalError);
    }
}


template<class Type>
typename Foam::CoeffField<Type>::scalarTypeField&
Foam::CoeffField<Type>::toScalar()
{
    if (!scalarCoeffPtr_)
    {
        // Debug check: demotion
        if (linearCoeffPtr_ || squareCoeffPtr_)
        {
            FatalErrorIn
            (
                "CoeffField<Type>::scalarTypeField& "
                "CoeffField<Type>::toScalar()"
            )   << "Detected demotion to scalar.  Probably an error"
                << abort(FatalError);
        }

        scalarCoeffPtr_ =
            new scalarTypeField(size(), pTraits<scalarType>::zero);
    }

    return *scalarCoeffPtr_;
}


template<class Type>
typename Foam::CoeffField<Type>::linearTypeField&
Foam::CoeffField<Type>::toLinear()
{
    if (!linearCoeffPtr_)
    {
        // Debug check: demotion
        if (squareCoeffPtr_)
        {
            FatalErrorIn
            (
                "CoeffField<Type>::linearTypeField& "
                "CoeffField<Type>::toLinear()"
            )   << "Detected demotion to linear.  Probably an error"
                << abort(FatalError);
        }

        linearCoeffPtr_ =
            new linearTypeField(size(), pTraits<linearType>::zero);

        // If scalar is active, promote to linear
        if (scalarCoeffPtr_)
        {
            expandScalar(*linearCoeffPtr_, *scalarCoeffPtr_);
            deleteDemandDrivenData(scalarCoeffPtr_);
        }
    }

    return *linearCoeffPtr_;
}


template<class Type>
typename Foam::CoeffField<Type>::squareTypeField&
Foam::CoeffField<Type>::toSquare()
{
    if (!squareCoeffPtr_)
    {
        squareCoeffPtr_ =
            new squareTypeField(this->size(), pTraits<squareType>::zero);

        // If scalar is active, promote to square
        if (scalarCoeffPtr_)
        {
            expandScalar(*squareCoeffPtr_, *scalarCoeffPtr_);
            deleteDemandDrivenData(scalarCoeffPtr_);
        }

        // If linear is active, promote to square
        if (linearCoeffPtr_)
        {
            expandLinear(*squareCoeffPtr_, *linearCoeffPtr_);
            deleteDemandDrivenData(linearCoeffPtr_);
        }
    }

    return *squareCoeffPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CoeffField<Type>::CoeffField(const label size)
:
    scalarCoeffPtr_(NULL),
    linearCoeffPtr_(NULL),
    squareCoeffPtr_(NULL),
    size_(size)
{}


template<class Type>
Foam::CoeffField<Type>::CoeffField(const CoeffField<Type>& f)
:
    refCount(),
    scalarCoeffPtr_(NULL),
    linearCoeffPtr_(NULL),
    squareCoeffPtr_(NULL),
    size_(f.size())
{
    if (f.scalarCoeffPtr_)
    {
        scalarCoeffPtr_ = new scalarTypeField(*(f.scalarCoeffPtr_));
    }
    else if (f.linearCoeffPtr_)
    {
        linearCoeffPtr_ = new linearTypeField(*(f.linearCoeffPtr_));
    }
    else if (f.squareCoeffPtr_)
    {
        squareCoeffPtr_ = new squareTypeField(*(f.squareCoeffPtr_));
    }
}


template<class Type>
Foam::CoeffField<Type>::CoeffField(Istream& is)
:
    scalarCoeffPtr_(NULL),
    linearCoeffPtr_(NULL),
    squareCoeffPtr_(NULL),
    size_(0)
{
    // Read keyword and pick up allocated field
    word key(is);

    if
    (
        key
     == blockCoeffBase::activeLevelNames_[blockCoeffBase::UNALLOCATED]
    )
    {
        size_ = readLabel(is);
    }
    else if
    (
        key
     == blockCoeffBase::activeLevelNames_[blockCoeffBase::SCALAR]
    )
    {
        scalarCoeffPtr_ = new scalarTypeField(is);
        size_ = scalarCoeffPtr_->size();
    }
    else if
    (
        key
     == blockCoeffBase::activeLevelNames_[blockCoeffBase::LINEAR]
    )
    {
        linearCoeffPtr_ = new linearTypeField(is);
        size_ = linearCoeffPtr_->size();
    }
    else if
    (
        key
     == blockCoeffBase::activeLevelNames_[blockCoeffBase::SQUARE]
    )
    {
        squareCoeffPtr_ = new squareTypeField(is);
        size_ = squareCoeffPtr_->size();
    }
    else
    {
        FatalIOErrorIn
        (
            "CoeffField<Type>::CoeffField(Istream& is)",
            is
        )   << "invalid keyword while reading: " << key
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::tmp<Foam::CoeffField<Type> > Foam::CoeffField<Type>::clone() const
{
    return tmp<CoeffField<Type> >(new CoeffField<Type>(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::CoeffField<Type>::~CoeffField()
{
    this->clear();
}


template<class Type>
void Foam::CoeffField<Type>::clear()
{
    deleteDemandDrivenData(scalarCoeffPtr_);
    deleteDemandDrivenData(linearCoeffPtr_);
    deleteDemandDrivenData(squareCoeffPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Foam::label Foam::CoeffField<Type>::size() const
{
    return size_;
}


template<class Type>
void Foam::CoeffField<Type>::negate()
{
    if (scalarCoeffPtr_)
    {
        scalarCoeffPtr_->negate();
    }
    else if (linearCoeffPtr_)
    {
        linearCoeffPtr_->negate();
    }
    else if (squareCoeffPtr_)
    {
        squareCoeffPtr_->negate();
    }
}


template<class Type>
Foam::tmp<Foam::CoeffField<Type> > Foam::CoeffField<Type>::transpose() const
{
    tmp<CoeffField<Type> > tt(new CoeffField<Type>(this->size()));
    CoeffField<Type>& t = tt();

    if (scalarCoeffPtr_)
    {
        t.toScalar() = *scalarCoeffPtr_;
    }
    else if (linearCoeffPtr_)
    {
        t.toLinear() = *linearCoeffPtr_;
    }
    else if (squareCoeffPtr_)
    {
        t.toSquare() = this->asSquare().T();
    }
    else
    {
        // Not allocated - do nothing
    }


    return tt;
}


template<class Type>
Foam::blockCoeffBase::activeLevel
Foam::CoeffField<Type>::activeType() const
{
    if (scalarCoeffPtr_)
    {
        return blockCoeffBase::SCALAR;
    }
    else if (linearCoeffPtr_)
    {
        return blockCoeffBase::LINEAR;
    }
    else if (squareCoeffPtr_)
    {
        return blockCoeffBase::SQUARE;
    }
    else
    {
        return blockCoeffBase::UNALLOCATED;
    }
}


template<class Type>
void Foam::CoeffField<Type>::checkActive() const
{
    label nActive = 0;

    if (scalarCoeffPtr_) nActive++;
    if (linearCoeffPtr_) nActive++;
    if (squareCoeffPtr_) nActive++;

    if (nActive > 1)
    {
        FatalErrorIn("void Foam::CoeffField<Type>::checkActive() const")
            << "Activation/deactivation error.  nActive = " << nActive
            << abort(FatalError);
    }
}


template<class Type>
const typename Foam::CoeffField<Type>::scalarTypeField&
Foam::CoeffField<Type>::asScalar() const
{
    if (!scalarCoeffPtr_)
    {
        FatalErrorIn
        (
            "CoeffField<Type>::scalarTypeField& "
            "CoeffField<Type>::asScalar()"
        )   << "Requested scalar but active type is: "
            << blockCoeffBase::activeLevelNames_[this->activeType()]
            << ".  This is not allowed."
            << abort(FatalError);
    }

    return *scalarCoeffPtr_;
}


template<class Type>
const typename Foam::CoeffField<Type>::linearTypeField&
Foam::CoeffField<Type>::asLinear() const
{
    if (!linearCoeffPtr_)
    {
        FatalErrorIn
        (
            "CoeffField<Type>::linearTypeField& "
            "CoeffField<Type>::asLinear()"
        )   << "Requested linear but active type is: "
            << blockCoeffBase::activeLevelNames_[this->activeType()]
            << ".  This is not allowed."
            << abort(FatalError);
    }

    return *linearCoeffPtr_;
}


template<class Type>
const typename Foam::CoeffField<Type>::squareTypeField&
Foam::CoeffField<Type>::asSquare() const
{
    if (!squareCoeffPtr_)
    {
        FatalErrorIn
        (
            "CoeffField<Type>::squareTypeField& "
            "CoeffField<Type>::asSquare()"
        )   << "Requested square but active type is: "
            << blockCoeffBase::activeLevelNames_[this->activeType()]
            << ".  This is not allowed."
            << abort(FatalError);
    }

    return *squareCoeffPtr_;
}


template<class Type>
typename Foam::CoeffField<Type>::scalarTypeField&
Foam::CoeffField<Type>::asScalar()
{
    if (linearCoeffPtr_ || squareCoeffPtr_)
    {
        FatalErrorIn
        (
            "CoeffField<Type>::scalarTypeField& "
            "CoeffField<Type>::asScalar()"
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
typename Foam::CoeffField<Type>::linearTypeField&
Foam::CoeffField<Type>::asLinear()
{
    if (squareCoeffPtr_)
    {
        FatalErrorIn
        (
            "CoeffField<Type>::linearTypeField& "
            "CoeffField<Type>::asLinear()"
        )   << "Requested linear but active type is: "
            << blockCoeffBase::activeLevelNames_[this->activeType()]
            << ".  This is not allowed."
            << abort(FatalError);
    }

    if (!linearCoeffPtr_)
    {
        return this->toLinear();
    }

    return *linearCoeffPtr_;
}


template<class Type>
typename Foam::CoeffField<Type>::squareTypeField&
Foam::CoeffField<Type>::asSquare()
{
    if (!squareCoeffPtr_)
    {
        return this->toSquare();
    }

    return *squareCoeffPtr_;
}


template<class Type>
Foam::tmp<typename Foam::CoeffField<Type>::scalarTypeField>
Foam::CoeffField<Type>::component(const direction dir) const
{
    if (scalarCoeffPtr_)
    {
        return *scalarCoeffPtr_;
    }
    else if (linearCoeffPtr_)
    {
        return linearCoeffPtr_->component(dir);
    }
    else if (squareCoeffPtr_)
    {
        linearTypeField lc(this->size());
        contractLinear(lc, *squareCoeffPtr_);

        return lc.component(dir);
    }
    else
    {
        FatalErrorIn
        (
            "tmp<CoeffField<Type>::scalarTypeField>"
            "CoeffField<Type>::component(const direction dir) const"
        )   << "Field not allocated."
            << abort(FatalError);
    }

    // Dummy return to keep compiler happy
    return *scalarCoeffPtr_;
}


template<class Type>
Foam::BlockCoeff<Type>
Foam::CoeffField<Type>::getCoeff(const label index) const
{
    BlockCoeff<Type> result;

    if (scalarCoeffPtr_)
    {
        result.asScalar() = (*scalarCoeffPtr_)[index];
    }
    else if (linearCoeffPtr_)
    {
        result.asLinear() = (*linearCoeffPtr_)[index];
    }
    else if (squareCoeffPtr_)
    {
        result.asSquare() = (*squareCoeffPtr_)[index];
    }

    return result;
}


template<class Type>
void Foam::CoeffField<Type>::setCoeff
(
    const label index,
    const BlockCoeff<Type>& coeff
)
{
    BlockCoeff<Type> result;

    if (coeff.activeType() == blockCoeffBase::SCALAR)
    {
        (*scalarCoeffPtr_)[index] = result.asScalar();
    }
    else if (coeff.activeType() == blockCoeffBase::LINEAR)
    {
        (*linearCoeffPtr_)[index] = result.asLinear();
    }
    else if (coeff.activeType() == blockCoeffBase::SQUARE)
    {
        (*squareCoeffPtr_)[index] = result.asSquare();
    }
}


template<class Type>
void Foam::CoeffField<Type>::getSubset
(
    CoeffField<Type>& f,
    const label start,
    const label size
) const
{
    // Check sizes
    if (f.size() != size)
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "void Foam::CoeffField<Type>::getSubset\n"
            "(\n"
            "    CoeffField<Type>& f,\n"
            "    const label start,\n"
            "    const label size\n"
            ") const"
        )   << "Incompatible sizes: " << f.size() << " and " << size
            << abort(FatalError);
    }

    if (scalarCoeffPtr_)
    {
        scalarTypeField& ff = f.asScalar();

        const scalarTypeField& localF = (*scalarCoeffPtr_);

        forAll (ff, ffI)
        {
            ff[ffI] = localF[start + ffI];
        }
    }
    else if (linearCoeffPtr_)
    {
        linearTypeField& ff = f.asLinear();

        const linearTypeField& localF = (*linearCoeffPtr_);

        forAll (ff, ffI)
        {
            ff[ffI] = localF[start + ffI];
        }
    }
    else if (squareCoeffPtr_)
    {
        squareTypeField& ff = f.asSquare();

        const squareTypeField& localF = (*squareCoeffPtr_);

        forAll (ff, ffI)
        {
            ff[ffI] = localF[start + ffI];
        }
    }
}


template<class Type>
void Foam::CoeffField<Type>::getSubset
(
    CoeffField<Type>& f,
    const labelList& addr
) const
{
    // Check sizes
    if (f.size() != addr.size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "void Foam::CoeffField<Type>::getSubset\n"
            "(\n"
            "    CoeffField<Type>& f,\n"
            "    const labelList addr\n"
            ") const"
        )   << "Incompatible sizes: " << f.size() << " and " << addr.size()
            << abort(FatalError);
    }

    if (scalarCoeffPtr_)
    {
        scalarTypeField& ff = f.asScalar();

        const scalarTypeField& localF = (*scalarCoeffPtr_);

        forAll (ff, ffI)
        {
            ff[ffI] = localF[addr[ffI]];
        }
    }
    else if (linearCoeffPtr_)
    {
        linearTypeField& ff = f.asLinear();

        const linearTypeField& localF = (*linearCoeffPtr_);

        forAll (ff, ffI)
        {
            ff[ffI] = localF[addr[ffI]];
        }
    }
    else if (squareCoeffPtr_)
    {
        squareTypeField& ff = f.asSquare();

        const squareTypeField& localF = (*squareCoeffPtr_);

        forAll (ff, ffI)
        {
            ff[ffI] = localF[addr[ffI]];
        }
    }
}


template<class Type>
void Foam::CoeffField<Type>::setSubset
(
    const CoeffField<Type>& f,
    const label start,
    const label size
)
{
    // Check sizes
    if (f.size() != size)
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "void Foam::CoeffField<Type>::setSubset\n"
            "(\n"
            "     const CoeffField<Type>& f,\n"
            "    const label start,\n"
            "    const label size\n"
            ")"
        )   << "Incompatible sizes: " << f.size() << " and " << size
            << abort(FatalError);
    }

    if (f.activeType() == blockCoeffBase::SCALAR)
    {
        const scalarTypeField& ff = f.asScalar();

        scalarTypeField& localF = this->asScalar();

        forAll (ff, ffI)
        {
            localF[start + ffI] = ff[ffI];
        }
    }
    else if (f.activeType() == blockCoeffBase::LINEAR)
    {
        const linearTypeField& ff = f.asLinear();

        linearTypeField& localF = this->asLinear();

        forAll (ff, ffI)
        {
            localF[start + ffI] = ff[ffI];
        }
    }
    else if (f.activeType() == blockCoeffBase::SQUARE)
    {
        const squareTypeField& ff = f.asSquare();

        squareTypeField& localF = this->asSquare();

        forAll (ff, ffI)
        {
            localF[start + ffI] = ff[ffI];
        }
    }
}


template<class Type>
void Foam::CoeffField<Type>::setSubset
(
    const CoeffField<Type>& f,
    const labelList& addr
)
{
    // Check sizes
    if (f.size() != addr.size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "void Foam::CoeffField<Type>::setSubset\n"
            "(\n"
            "    const CoeffField<Type>& f,\n"
            "    const labelList addr\n"
            ")"
        )   << "Incompatible sizes: " << f.size() << " and " << addr.size()
            << abort(FatalError);
    }

    if (f.activeType() == blockCoeffBase::SCALAR)
    {
        const scalarTypeField& ff = f.asScalar();

        scalarTypeField& localF = this->asScalar();

        forAll (ff, ffI)
        {
            localF[addr[ffI]] = ff[ffI];
        }
    }
    else if (f.activeType() == blockCoeffBase::LINEAR)
    {
        const linearTypeField& ff = f.asLinear();

        linearTypeField& localF = this->asLinear();

        forAll (ff, ffI)
        {
            localF[addr[ffI]] = ff[ffI];
        }
    }
    else if (f.activeType() == blockCoeffBase::SQUARE)
    {
        const squareTypeField& ff = f.asSquare();

        squareTypeField& localF = this->asSquare();

        forAll (ff, ffI)
        {
            localF[addr[ffI]] = ff[ffI];
        }
    }
}


template<class Type>
void Foam::CoeffField<Type>::zeroOutSubset
(
    const label start,
    const label size
)
{
    if (scalarCoeffPtr_)
    {
        scalarTypeField& localF = this->asScalar();

        for (label ffI = 0; ffI < size; ffI++)
        {
            localF[start + ffI] = pTraits<scalarType>::zero;
        }
    }
    else if (linearCoeffPtr_)
    {
        linearTypeField& localF = this->asLinear();

        for (label ffI = 0; ffI < size; ffI++)
        {
            localF[start + ffI] = pTraits<linearType>::zero;
        }
    }
    else if (squareCoeffPtr_)
    {
        squareTypeField& localF = this->asSquare();

        for (label ffI = 0; ffI < size; ffI++)
        {
            localF[start + ffI] = pTraits<squareType>::zero;
        }
    }
}


template<class Type>
void Foam::CoeffField<Type>::zeroOutSubset
(
    const labelList& addr
)
{
    if (scalarCoeffPtr_)
    {
        scalarTypeField& localF = this->asScalar();

        forAll (addr, ffI)
        {
            localF[addr[ffI]] = pTraits<scalarType>::zero;
        }
    }
    else if (linearCoeffPtr_)
    {
        linearTypeField& localF = this->asLinear();

        forAll (addr, ffI)
        {
            localF[addr[ffI]] = pTraits<linearType>::zero;
        }
    }
    else if (squareCoeffPtr_)
    {
        squareTypeField& localF = this->asSquare();

        forAll (addr, ffI)
        {
            localF[addr[ffI]] = pTraits<squareType>::zero;
        }
    }
}


template<class Type>
void Foam::CoeffField<Type>::addSubset
(
    const CoeffField<Type>& f,
    const labelList& addr
)
{
    // Check sizes
    if (f.size() != addr.size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "void Foam::CoeffField<Type>::addSubset\n"
            "(\n"
            "    const CoeffField<Type>& f,\n"
            "    const labelList addr\n"
            ")"
        )   << "Incompatible sizes: " << f.size() << " and " << addr.size()
            << abort(FatalError);
    }

    if (this->activeType() == blockCoeffBase::SQUARE)
    {
        if (f.activeType() == blockCoeffBase::SQUARE)
        {
            squareTypeField& localF = this->asSquare();

            const squareTypeField& ff = f.asSquare();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] += ff[ffI];
            }
        }
        else if (f.activeType() == blockCoeffBase::LINEAR)
        {
            squareTypeField& localF = this->asSquare();

            squareTypeField ff(f.size());
            expandLinear(ff, f.asLinear());

            forAll (ff, ffI)
            {
                localF[addr[ffI]] += ff[ffI];
            }
        }
        else if (f.activeType() == blockCoeffBase::SCALAR)
        {
            squareTypeField& localF = this->asSquare();

            squareTypeField ff(f.size());
            expandScalar(ff, f.asScalar());

            forAll (ff, ffI)
            {
                localF[addr[ffI]] += ff[ffI];
            }
        }
    }
    else if (this->activeType() == blockCoeffBase::LINEAR)
    {
        if (f.activeType() == blockCoeffBase::SQUARE)
        {
            squareTypeField& localF = this->asSquare();

            const squareTypeField& ff = f.asSquare();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] += ff[ffI];
            }
        }
        if (f.activeType() == blockCoeffBase::LINEAR)
        {
            linearTypeField& localF = this->asLinear();

            const linearTypeField& ff = f.asLinear();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] += ff[ffI];
            }
        }
        else if (f.activeType() == blockCoeffBase::SCALAR)
        {
            linearTypeField& localF = this->asLinear();

            linearTypeField ff(f.size());
            expandScalar(ff, f.asScalar());

            forAll (ff, ffI)
            {
                localF[addr[ffI]] += ff[ffI];
            }
        }
    }
    else if (this->activeType() == blockCoeffBase::SCALAR)
    {
        if (f.activeType() == blockCoeffBase::SQUARE)
        {
            squareTypeField& localF = this->asSquare();

            const squareTypeField& ff = f.asSquare();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] += ff[ffI];
            }
        }
        if (f.activeType() == blockCoeffBase::LINEAR)
        {
            linearTypeField& localF = this->asLinear();

            const linearTypeField& ff = f.asLinear();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] += ff[ffI];
            }
        }
        else if (f.activeType() == blockCoeffBase::SCALAR)
        {
            const scalarTypeField& ff = f.asScalar();

            scalarTypeField& localF = this->asScalar();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] += ff[ffI];
            }
        }
    }
    else
    {
        if (f.activeType() == blockCoeffBase::SQUARE)
        {
            squareTypeField& localF = this->toSquare();

            const squareTypeField& ff = f.asSquare();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] += ff[ffI];
            }
        }
        if (f.activeType() == blockCoeffBase::LINEAR)
        {
            linearTypeField& localF = this->toLinear();

            const linearTypeField& ff = f.asLinear();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] += ff[ffI];
            }
        }
        else if (f.activeType() == blockCoeffBase::SCALAR)
        {
            scalarTypeField& localF = this->toScalar();

            const scalarTypeField& ff = f.asScalar();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] += ff[ffI];
            }
        }
//         FatalErrorIn
//         (
//             "template<class Type>\n"
//             "void Foam::CoeffField<Type>::addSubset\n"
//             "(\n"
//             "    const CoeffField<Type>& f,\n"
//             "    const labelList addr\n"
//             ")"
//         )   << "Incompatible combination of types"
//             << abort(FatalError);
    }
}


template<class Type>
void Foam::CoeffField<Type>::subtractSubset
(
    const CoeffField<Type>& f,
    const labelList& addr
)
{
    // Check sizes
    if (f.size() != addr.size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "void Foam::CoeffField<Type>::subtractSubset\n"
            "(\n"
            "    const CoeffField<Type>& f,\n"
            "    const labelList addr\n"
            ")"
        )   << "Incompatible sizes: " << f.size() << " and " << addr.size()
            << abort(FatalError);
    }

    if (this->activeType() == blockCoeffBase::SQUARE)
    {
        if (f.activeType() == blockCoeffBase::SQUARE)
        {
            squareTypeField& localF = this->asSquare();

            const squareTypeField& ff = f.asSquare();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] -= ff[ffI];
            }
        }
        else if (f.activeType() == blockCoeffBase::LINEAR)
        {
            squareTypeField& localF = this->asSquare();

            squareTypeField ff(f.size());
            expandLinear(ff, f.asLinear());

            forAll (ff, ffI)
            {
                localF[addr[ffI]] -= ff[ffI];
            }
        }
        else if (f.activeType() == blockCoeffBase::SCALAR)
        {
            squareTypeField& localF = this->asSquare();

            squareTypeField ff(f.size());
            expandScalar(ff, f.asScalar());

            forAll (ff, ffI)
            {
                localF[addr[ffI]] -= ff[ffI];
            }
        }
    }
    else if (this->activeType() == blockCoeffBase::LINEAR)
    {
        if (f.activeType() == blockCoeffBase::SQUARE)
        {
            squareTypeField& localF = this->asSquare();

            const squareTypeField& ff = f.asSquare();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] -= ff[ffI];
            }
        }
        if (f.activeType() == blockCoeffBase::LINEAR)
        {
            linearTypeField& localF = this->asLinear();

            const linearTypeField& ff = f.asLinear();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] -= ff[ffI];
            }
        }
        else if (f.activeType() == blockCoeffBase::SCALAR)
        {
            linearTypeField& localF = this->asLinear();

            linearTypeField ff(f.size());
            expandScalar(ff, f.asScalar());

            forAll (ff, ffI)
            {
                localF[addr[ffI]] -= ff[ffI];
            }
        }
    }
    else if (this->activeType() == blockCoeffBase::SCALAR)
    {
        if (f.activeType() == blockCoeffBase::SQUARE)
        {
            squareTypeField& localF = this->asSquare();

            const squareTypeField& ff = f.asSquare();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] -= ff[ffI];
            }
        }
        if (f.activeType() == blockCoeffBase::LINEAR)
        {
            linearTypeField& localF = this->asLinear();

            const linearTypeField& ff = f.asLinear();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] -= ff[ffI];
            }
        }
        else if (f.activeType() == blockCoeffBase::SCALAR)
        {
            const scalarTypeField& ff = f.asScalar();

            scalarTypeField& localF = this->asScalar();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] -= ff[ffI];
            }
        }
    }
    else
    {
        if (f.activeType() == blockCoeffBase::SQUARE)
        {
            squareTypeField& localF = this->toSquare();

            const squareTypeField& ff = f.asSquare();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] -= ff[ffI];
            }
        }
        if (f.activeType() == blockCoeffBase::LINEAR)
        {
            linearTypeField& localF = this->toLinear();

            const linearTypeField& ff = f.asLinear();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] -= ff[ffI];
            }
        }
        else if (f.activeType() == blockCoeffBase::SCALAR)
        {
            scalarTypeField& localF = this->toScalar();

            const scalarTypeField& ff = f.asScalar();

            forAll (ff, ffI)
            {
                localF[addr[ffI]] -= ff[ffI];
            }
        }
//         FatalErrorIn
//         (
//             "template<class Type>\n"
//             "void Foam::CoeffField<Type>::subtractSubset\n"
//             "(\n"
//             "    const CoeffField<Type>& f,\n"
//             "    const labelList addr\n"
//             ")"
//         )   << "Incompatible combination of types"
//             << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::CoeffField<Type>::operator=(const CoeffField<Type>& f)
{
    if (this == &f)
    {
        FatalErrorIn("CoeffField<Type>::operator=(const CoeffField<Type>&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    // Check field sizes
    if (f.size() != this->size())
    {
        FatalErrorIn
        (
            "void CoeffField<Type>::operator=(const CoeffField<Type>& f)"
        )   << "Incorrect field size: " << f.size()
            << " local size: " << size()
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
    else if (f.squareCoeffPtr_)
    {
        this->toSquare() = *(f.squareCoeffPtr_);
    }
    else
    {
        // Not allocated - do nothing
    }
}


template<class Type>
void Foam::CoeffField<Type>::operator=(const tmp<CoeffField<Type> >& tf)
{
    if (this == &(tf()))
    {
        FatalErrorIn("CoeffField<Type>::operator=(const tmp<CoeffField>&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    operator=(tf());
    tf.clear();
}


template<class Type>
Foam::tmp<Foam::CoeffField<Type> > Foam::CoeffField<Type>::operator-()
{
    tmp<CoeffField<Type> > tf(new CoeffField<Type>(*this));
    tf->negate();
    return tf;
}


#define COMPUTED_BASE_ASSIGNMENT(op,defaultOp)                                \
                                                                              \
template<class Type>                                                          \
void Foam::CoeffField<Type>::operator op(const CoeffField<Type>& f)           \
{                                                                             \
    if (f.size() != this->size())                                             \
    {                                                                         \
        FatalErrorIn                                                          \
        (                                                                     \
            "void CoeffField<tensor>::operator "                              \
            "op(const CoeffField<tensor>& f)"                                 \
        )   << "Incorrect field size: " << f.size()                           \
            << " local size: " << size()                                      \
            << abort(FatalError);                                             \
    }                                                                         \
                                                                              \
    if (this->activeType() == blockCoeffBase::SQUARE)                         \
    {                                                                         \
        if (f.activeType() == blockCoeffBase::SQUARE)                         \
        {                                                                     \
            this->asSquare() op f.asSquare();                                 \
        }                                                                     \
        else if (f.activeType() == blockCoeffBase::LINEAR)                    \
        {                                                                     \
            squareTypeField ff(f.size());                                     \
            expandLinear(ff, f.asLinear());                                   \
                                                                              \
            this->asSquare() op ff;                                           \
        }                                                                     \
        else if (f.activeType() == blockCoeffBase::SCALAR)                    \
        {                                                                     \
            squareTypeField ff(f.size());                                     \
            expandScalar(ff, f.asScalar());                                   \
                                                                              \
            this->asSquare() op ff;                                           \
        }                                                                     \
    }                                                                         \
    else if (this->activeType() == blockCoeffBase::LINEAR)                    \
    {                                                                         \
        if (f.activeType() == blockCoeffBase::SQUARE)                         \
        {                                                                     \
            this->asSquare() op f.asSquare();                                 \
        }                                                                     \
        if (f.activeType() == blockCoeffBase::LINEAR)                         \
        {                                                                     \
            this->asLinear() op f.asLinear();                                 \
        }                                                                     \
        else if (f.activeType() == blockCoeffBase::SCALAR)                    \
        {                                                                     \
            linearTypeField ff(f.size());                                     \
            expandScalar(ff, f.asScalar());                                   \
                                                                              \
            this->asLinear() op ff;                                           \
        }                                                                     \
    }                                                                         \
    else if (this->activeType() == blockCoeffBase::SCALAR)                    \
    {                                                                         \
        if (f.activeType() == blockCoeffBase::SQUARE)                         \
        {                                                                     \
            this->asSquare() op f.asSquare();                                 \
        }                                                                     \
        if (f.activeType() == blockCoeffBase::LINEAR)                         \
        {                                                                     \
            this->asLinear() op f.asLinear();                                 \
        }                                                                     \
        else if (f.activeType() == blockCoeffBase::SCALAR)                    \
        {                                                                     \
            this->asScalar() op f.asScalar();                                 \
        }                                                                     \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        if (f.activeType() == blockCoeffBase::SQUARE)                         \
        {                                                                     \
            (*this) defaultOp f.asSquare();                                   \
        }                                                                     \
        if (f.activeType() == blockCoeffBase::LINEAR)                         \
        {                                                                     \
            (*this) defaultOp f.asLinear();                                   \
        }                                                                     \
        else if (f.activeType() == blockCoeffBase::SCALAR)                    \
        {                                                                     \
            (*this) defaultOp f.asScalar();                                   \
        }                                                                     \
    }                                                                         \
}                                                                             \
                                                                              \
template<class Type>                                                          \
void Foam::CoeffField<Type>::operator op(const tmp<CoeffField<Type> >& tf)    \
{                                                                             \
    operator op(tf());                                                        \
    tf.clear();                                                               \
}


#define COMPUTED_ARG_ASSIGNMENT(op)                                           \
                                                                              \
template<class Type>                                                          \
void Foam::CoeffField<Type>::operator op(const scalarTypeField& f)            \
{                                                                             \
    checkSize(f);                                                             \
                                                                              \
    const blockCoeffBase::activeLevel al = this->activeType();                \
                                                                              \
    if (al == blockCoeffBase::UNALLOCATED || al == blockCoeffBase::SCALAR)    \
    {                                                                         \
        this->toScalar() op f;                                                \
    }                                                                         \
    else if (al == blockCoeffBase::LINEAR)                                    \
    {                                                                         \
        this->toLinear() op f*pTraits<linearType>::one;                       \
    }                                                                         \
    else if (al == blockCoeffBase::SQUARE)                                    \
    {                                                                         \
        squareTypeField stf(f.size());                                        \
        expandScalar(stf, f);                                                 \
        this->toSquare() op stf;                                              \
    }                                                                         \
}                                                                             \
                                                                              \
template<class Type>                                                          \
void Foam::CoeffField<Type>::operator op(const tmp<scalarTypeField>& tf)      \
{                                                                             \
    operator op(tf());                                                        \
    tf.clear();                                                               \
}                                                                             \
                                                                              \
                                                                              \
template<class Type>                                                          \
void Foam::CoeffField<Type>::operator op(const linearTypeField& f)            \
{                                                                             \
    checkSize(f);                                                             \
                                                                              \
    const blockCoeffBase::activeLevel al = this->activeType();                \
                                                                              \
    if                                                                        \
    (                                                                         \
        al == blockCoeffBase::UNALLOCATED                                     \
     || al == blockCoeffBase::SCALAR                                          \
     || al == blockCoeffBase::LINEAR                                          \
    )                                                                         \
    {                                                                         \
        this->toLinear() op f;                                                \
    }                                                                         \
    else if (al == blockCoeffBase::SQUARE)                                    \
    {                                                                         \
        squareTypeField stf(f.size());                                        \
        expandLinear(stf, f);                                                 \
        this->toSquare() op stf;                                              \
    }                                                                         \
}                                                                             \
                                                                              \
template<class Type>                                                          \
void Foam::CoeffField<Type>::operator op(const tmp<linearTypeField>& tf)      \
{                                                                             \
    operator op(tf());                                                        \
    tf.clear();                                                               \
}                                                                             \
                                                                              \
                                                                              \
template<class Type>                                                          \
void Foam::CoeffField<Type>::operator op(const squareTypeField& f)            \
{                                                                             \
    checkSize(f);                                                             \
    this->toSquare() op f;                                                    \
}                                                                             \
                                                                              \
template<class Type>                                                          \
void Foam::CoeffField<Type>::operator op(const tmp<squareTypeField>& tf)      \
{                                                                             \
    operator op(tf());                                                        \
    tf.clear();                                                               \
}

#define COMPUTED_BASE_OPERATOR(TYPE, op)                                      \
                                                                              \
template<class Type>                                                          \
void Foam::CoeffField<Type>::operator op(const TYPE& t)                       \
{                                                                             \
    if (scalarCoeffPtr_)                                                      \
    {                                                                         \
        *(scalarCoeffPtr_) op t;                                              \
    }                                                                         \
    else if (linearCoeffPtr_)                                                 \
    {                                                                         \
        *(linearCoeffPtr_) op t;                                              \
    }                                                                         \
    else if (squareCoeffPtr_)                                                 \
    {                                                                         \
        *(squareCoeffPtr_) op t;                                              \
    }                                                                         \
    else                                                                      \
    {                                                                         \
    }                                                                         \
}                                                                             \
                                                                              \
template<class Type>                                                          \
void Foam::CoeffField<Type>::operator op(const UList<TYPE>& tf)               \
{                                                                             \
    checkSize(tf);                                                            \
                                                                              \
    if (scalarCoeffPtr_)                                                      \
    {                                                                         \
        *(scalarCoeffPtr_) op tf;                                             \
    }                                                                         \
    else if (linearCoeffPtr_)                                                 \
    {                                                                         \
        *(linearCoeffPtr_) op tf;                                             \
    }                                                                         \
    else if (squareCoeffPtr_)                                                 \
    {                                                                         \
        *(squareCoeffPtr_) op tf;                                             \
    }                                                                         \
    else                                                                      \
    {                                                                         \
    }                                                                         \
}                                                                             \
                                                                              \
template<class Type>                                                          \
void Foam::CoeffField<Type>::operator op(const tmp<Field<TYPE> >& tf)         \
{                                                                             \
    operator op(tf());                                                        \
    tf.clear();                                                               \
}


#define COMPUTED_ASSIGNMENT(op, defaultOp)                                    \
COMPUTED_BASE_ASSIGNMENT(op, defaultOp)                                       \
COMPUTED_ARG_ASSIGNMENT(op)

// Remaining operator=
COMPUTED_ARG_ASSIGNMENT(=)

COMPUTED_ASSIGNMENT(+=, =)
COMPUTED_ASSIGNMENT(-=, =-)

COMPUTED_BASE_OPERATOR(scalar, *=)
COMPUTED_BASE_OPERATOR(scalar, /=)

#undef COMPUTED_BASE_OPERATOR
#undef COMPUTED_BASE_ASSIGNMENT
#undef COMPUTED_ARG_ASSIGNMENT
#undef COMPUTED_ASSIGNMENT


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const CoeffField<Type>& f)
{
    // Write active type
    os << blockCoeffBase::activeLevelNames_[f.activeType()] << nl;

    if (f.activeType() == blockCoeffBase::SCALAR)
    {
        os << f.asScalar();
    }
    else if (f.activeType() == blockCoeffBase::LINEAR)
    {
        os << f.asLinear();
    }
    else if (f.activeType() == blockCoeffBase::SQUARE)
    {
        os << f.asSquare();
    }
    else
    {
        // Not allocated: write size
        os << f.size();
    }

    return os;
}


template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const tmp<CoeffField<Type> >& tf)
{
    os << tf();
    tf.clear();
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "CoeffFieldFunctions.C"


// ************************************************************************* //
