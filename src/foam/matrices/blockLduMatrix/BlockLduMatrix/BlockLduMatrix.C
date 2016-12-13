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

Description
    BlockLduMatrix is a general matrix class in which the coefficients are
    stored as three arrays, one for the upper triangle, one for the
    lower triangle and a third for the diagonal.  Addressing object must
    be supplied for the upper and lower triangles.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "BlockLduMatrix.H"
#include "IOstreams.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::debug::optimisationSwitch
Foam::BlockLduMatrix<Type>::fixFillIn
(
    "matrixConstraintFillIn",
    4
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockLduMatrix<Type>::BlockLduMatrix(const lduMesh& ldu)
:
    lduMesh_(ldu),
    diagPtr_(NULL),
    upperPtr_(NULL),
    lowerPtr_(NULL),
    interfaces_(ldu.interfaces().size()),
    coupleUpper_(ldu.lduAddr().nPatches()),
    coupleLower_(ldu.lduAddr().nPatches()),
    fixedEqns_(ldu.lduAddr().size()/fixFillIn())
{
    const lduAddressing& addr = ldu.lduAddr();

    forAll (coupleUpper_, i)
    {
        coupleUpper_.set(i, new CoeffField<Type>(addr.patchAddr(i).size()));
        coupleLower_.set(i, new CoeffField<Type>(addr.patchAddr(i).size()));
    }

}


template<class Type>
Foam::BlockLduMatrix<Type>::BlockLduMatrix(const BlockLduMatrix<Type>& A)
:
    refCount(),
    lduMesh_(A.lduMesh_),
    diagPtr_(NULL),
    upperPtr_(NULL),
    lowerPtr_(NULL),
    interfaces_(A.interfaces_),
    coupleUpper_(A.coupleUpper_),
    coupleLower_(A.coupleLower_),
    fixedEqns_(A.fixedEqns_)
{
    if (A.diagPtr_)
    {
        diagPtr_ = new TypeCoeffField(*(A.diagPtr_));
    }

    if (A.upperPtr_)
    {
        upperPtr_ = new TypeCoeffField(*(A.upperPtr_));
    }

    if (A.lowerPtr_)
    {
        lowerPtr_ = new TypeCoeffField(*(A.lowerPtr_));
    }
}


//HJ, problematic: memory management.
// Reconsider.  HJ, 7/Nov/2010
template<class Type>
Foam::BlockLduMatrix<Type>::BlockLduMatrix(BlockLduMatrix<Type>& A, bool reUse)
:
    refCount(),
    lduMesh_(A.lduMesh_),
    diagPtr_(NULL),
    upperPtr_(NULL),
    lowerPtr_(NULL),
    interfaces_(A.interfaces_, reUse),
    coupleUpper_(A.coupleUpper_, reUse),
    coupleLower_(A.coupleLower_, reUse),
    fixedEqns_(A.fixedEqns_)
{
    if (reUse)
    {
        if (A.lowerPtr_)
        {
            lowerPtr_ = A.lowerPtr_;
            A.lowerPtr_ = NULL;
        }

        if (A.diagPtr_)
        {
            diagPtr_ = A.diagPtr_;
            A.diagPtr_ = NULL;
        }

        if (A.upperPtr_)
        {
            upperPtr_ = A.upperPtr_;
            A.upperPtr_ = NULL;
        }
    }
    else
    {
        if (A.diagPtr_)
        {
            diagPtr_ = new TypeCoeffField(*(A.diagPtr_));
        }

        if (A.upperPtr_)
        {
            upperPtr_ = new TypeCoeffField(*(A.upperPtr_));
        }

        if (A.lowerPtr_)
        {
            lowerPtr_ = new TypeCoeffField(*(A.lowerPtr_));
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockLduMatrix<Type>::~BlockLduMatrix()
{
    deleteDemandDrivenData(diagPtr_);
    deleteDemandDrivenData(upperPtr_);
    deleteDemandDrivenData(lowerPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
typename Foam::BlockLduMatrix<Type>::TypeCoeffField&
Foam::BlockLduMatrix<Type>::diag()
{
    if (!diagPtr_)
    {
        diagPtr_ = new TypeCoeffField(lduAddr().size());
    }

    return *diagPtr_;
}


template<class Type>
const typename Foam::BlockLduMatrix<Type>::TypeCoeffField&
Foam::BlockLduMatrix<Type>::diag() const
{
    if (!diagPtr_)
    {
        FatalErrorIn
        (
            "const TypeCoeffField& BlockLduMatrix<Type>::diag() const"
        )   << "diagPtr_ unallocated"
            << abort(FatalError);
    }

    return *diagPtr_;
}


template<class Type>
typename Foam::BlockLduMatrix<Type>::TypeCoeffField&
Foam::BlockLduMatrix<Type>::upper()
{
    if (!upperPtr_)
    {
        upperPtr_ = new TypeCoeffField(lduAddr().lowerAddr().size());
    }

    return *upperPtr_;
}


template<class Type>
const typename Foam::BlockLduMatrix<Type>::TypeCoeffField&
Foam::BlockLduMatrix<Type>::upper() const
{
    if (!upperPtr_)
    {
        FatalErrorIn
        (
            "const TypeCoeffField& BlockLduMatrix<Type>::upper() const"
        )   << "upperPtr_ unallocated"
            << abort(FatalError);
    }

    return *upperPtr_;
}


template<class Type>
typename Foam::BlockLduMatrix<Type>::TypeCoeffField&
Foam::BlockLduMatrix<Type>::lower()
{
    if (!lowerPtr_)
    {
        if (upperPtr_)
        {
//             Info << "Manufacturing lower from upper transpose" << endl;
            lowerPtr_ = new TypeCoeffField(upperPtr_->transpose());
        }
        else
        {
            lowerPtr_ = new TypeCoeffField(lduAddr().lowerAddr().size());
        }
    }

    return *lowerPtr_;
}


template<class Type>
const typename Foam::BlockLduMatrix<Type>::TypeCoeffField&
Foam::BlockLduMatrix<Type>::lower() const
{
    if (!lowerPtr_)
    {
        FatalErrorIn
        (
            "const TypeCoeffField&  BlockLduMatrix<Type>::lower() const"
        )   << "lowerPtr_ unallocated"
            << abort(FatalError);
    }

    return *lowerPtr_;
}


template<class Type>
void Foam::BlockLduMatrix<Type>::clearInterfaces()
{
    forAll (interfaces_, i)
    {
        if (interfaces_.set(i))
        {
            delete interfaces_(i);
        }
    }
}



template<class Type>
bool Foam::BlockLduMatrix<Type>::empty() const
{
    return (!diagPtr_ && !lowerPtr_ && !upperPtr_);
}


template<class Type>
bool Foam::BlockLduMatrix<Type>::diagonal() const
{
    return (diagPtr_ && !lowerPtr_ && !upperPtr_);
}


template<class Type>
bool Foam::BlockLduMatrix<Type>::symmetric() const
{
    if (lowerPtr_ && !upperPtr_)
    {
        FatalErrorIn
        (
            "bool BlockLduMatrix<Type>::symmetric() const"
        )   << "Matrix assembly error: symmetric matrix but only lower "
            << "triangle is allocated.  This is not allowed."
            << abort(FatalError);
    }

    // Note:
    // It is possible that the block matrix is symmetric in sparseness
    // pattern, but asymmetric because the diagonal or upper coefficients are
    // square and asymmetric within the coefficients.
    // In such cases, the symmetric check in this function will falsely claim
    // the matrix is symmetries whereas in its flattened nature it is not
    // For the moment, symmetric check is used to see if the lower() is
    // allocated or not.
    // Please reconsider, especially related to matrix structure in
    // complex cases and in choice of linear equation solver for matrices
    // with asymmetric square coefficients
    // HJ and VV, 16/Mar/2016

    return (diagPtr_ && (!lowerPtr_ && upperPtr_));
}


template<class Type>
bool Foam::BlockLduMatrix<Type>::asymmetric() const
{
    return (diagPtr_ && lowerPtr_ && upperPtr_);
}


template<class Type>
bool Foam::BlockLduMatrix<Type>::componentCoupled() const
{
    // Return true if the matrix coefficient couple the components
    if (thereIsDiag())
    {
        if (diag().activeType() == blockCoeffBase::SQUARE)
        {
            return true;
        }
    }

    if (thereIsUpper())
    {
        if (upper().activeType() == blockCoeffBase::SQUARE)
        {
            return true;
        }
    }

    if (thereIsLower())
    {
        if (lower().activeType() == blockCoeffBase::SQUARE)
        {
            return true;
        }
    }

    return false;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const BlockLduMatrix<Type>& ldum)
{
    if (ldum.diagPtr_)
    {
        os.writeKeyword("diagonal")
            << *ldum.diagPtr_ <<  token::END_STATEMENT << nl;
    }
    else
    {
        // Dummy write for consistency
        os.writeKeyword("diagonal")
            << typename BlockLduMatrix<Type>::TypeCoeffField
               (
                   ldum.lduAddr().size()
               )
            << token::END_STATEMENT << nl;
    }

    if (ldum.upperPtr_)
    {
        os.writeKeyword("upper")
            << *ldum.upperPtr_
            << token::END_STATEMENT << nl;
    }
    else
    {
        // Dummy write for consistency
        os.writeKeyword("upper")
            << typename BlockLduMatrix<Type>::TypeCoeffField
               (
                   ldum.lduAddr().lowerAddr().size()
               )
            <<  token::END_STATEMENT << nl;
    }

    if (ldum.lowerPtr_)
    {
        os.writeKeyword("lower")
            << *ldum.lowerPtr_ <<  token::END_STATEMENT << nl;
    }
    else
    {
        // Dummy write for consistency
        os.writeKeyword("lower")
            << typename BlockLduMatrix<Type>::TypeCoeffField
               (
                   ldum.lduAddr().lowerAddr().size()
               )
            <<  token::END_STATEMENT << nl;
    }

    os.writeKeyword("coupleUpper")
        << ldum.coupleUpper_
        << token::END_STATEMENT << endl;

    os.writeKeyword("coupleLower")
        << ldum.coupleLower_
        << token::END_STATEMENT << endl;

    os.check("Ostream& operator<<(Ostream&, const BlockLduMatrix<Type>&");

    return os;
}


// ************************************************************************* //
