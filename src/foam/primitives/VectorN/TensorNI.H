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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template <class Cmpt, int length>
const char* const TensorN<Cmpt, length>::typeName =
    ("tensor" + name(length)).c_str();

template <class Cmpt, int length>
const TensorN<Cmpt, length> TensorN<Cmpt, length>::zero(0);

template <class Cmpt, int length>
const TensorN<Cmpt, length> TensorN<Cmpt, length>::one(1);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
template <class Cmpt, int length>
inline TensorN<Cmpt, length>::TensorN()
{}


//- Construct given VectorSpace
template <class Cmpt, int length>
inline TensorN<Cmpt, length>::TensorN
(
    const VectorSpace<TensorN<Cmpt, length>, Cmpt, length*length>& vs
)
:
    VectorSpace<TensorN<Cmpt, length>, Cmpt, length*length>(vs)
{}


//- Construct from component
template <class Cmpt, int length>
inline TensorN<Cmpt, length>::TensorN(const Cmpt& tx)
{
    VectorSpaceOps<TensorN<Cmpt, length>::nComponents,0>::eqOpS
    (
        *this,
        tx,
        eqOp<Cmpt>()
    );
}


//- Construct from Istream
template <class Cmpt, int length>
inline TensorN<Cmpt, length>::TensorN(Istream& is)
:
    VectorSpace<TensorN<Cmpt, length>, Cmpt, length*length>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return tensor diagonal
template <class Cmpt, int length>
inline DiagTensorN<Cmpt, length> TensorN<Cmpt, length>::diag() const
{
    DiagTensorN<Cmpt, length> dt;

    label diagI = 0;
    for (label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        dt[i] = this->v_[diagI];
        diagI += TensorN<Cmpt, length>::rowLength + 1;
    }

    return dt;
}


//- Return tensor transpose
template <class Cmpt, int length>
inline TensorN<Cmpt, length> TensorN<Cmpt, length>::T() const
{
    TensorN<Cmpt, length> transpose;

    label i = 0;
    for (label row = 0; row < TensorN<Cmpt, length>::rowLength; row++)
    {
        label j = row;
        for (label col = 0; col < TensorN<Cmpt, length>::rowLength; col++)
        {
            transpose.v_[i] = this->v_[j];
            i++;
            j += TensorN<Cmpt, length>::rowLength;
        }
    }

    return transpose;
}


//- Negative sum the vertical off-diagonal components
template <class Cmpt, int length>
inline TensorN<Cmpt, length> TensorN<Cmpt, length>::negSumDiag() const
{
    TensorN<Cmpt, length> negsumdiag;

    // Zero main diagonal
    label diagI = 0;
    for (label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        negsumdiag.v_[diagI] = 0.0;
        diagI += TensorN<Cmpt, length>::rowLength + 1;
    }

    label k = 0;
    for (label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        diagI = 0;
        for (label j = 0; j < TensorN<Cmpt, length>::rowLength; j++)
        {
            if (k != diagI)
            {
                negsumdiag.v_[k] = this->v_[k];
                negsumdiag.v_[diagI] -= this->v_[k];
            }
            k++;
            diagI += TensorN<Cmpt, length>::rowLength + 1;
        }
    }

    return negsumdiag;
}


//- Assign to a SphericalTensorN
template <class Cmpt, int length>
inline void
TensorN<Cmpt, length>::operator=(const SphericalTensorN<Cmpt, length>& st)
{
    label diag = 0;
    for (label i = 0; i < TensorN<Cmpt, length>::nComponents; i++)
    {
        if (i == diag)
        {
            this->v_[i] = st[0];
            diag += TensorN<Cmpt, length>::rowLength + 1;
        }
        else
        {
            this->v_[i] = pTraits<Cmpt>::zero;
        }
    }
}


//- Assign to a DiagTensorN
template <class Cmpt, int length>
inline void
TensorN<Cmpt, length>::operator=(const DiagTensorN<Cmpt, length>& dt)
{
    label diag = 0;
    label k = 0;
    for (label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        for (label j = 0; j < TensorN<Cmpt, length>::rowLength; j++)
        {
            if (j == diag)
            {
                this->v_[k] = dt[i];
            }
            else
            {
                this->v_[k] = pTraits<Cmpt>::zero;
            }
            k++;
        }
        diag++;
    }
}


//- Transform the tensor
//- The components are assumed to be individual scalars
//- i.e. transform has no effect
template<class Cmpt, int length>
inline TensorN<Cmpt, length> transform
(
    const tensor& tt,
    const TensorN<Cmpt, length>& v
)
{
    return v;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template <class Cmpt, int length>
inline direction TensorN<Cmpt, length>::cmpt
(
    const direction i,
    const direction j
)
{
    if
    (
        i > TensorN<Cmpt, length>::rowLength - 1
     || j > TensorN<Cmpt, length>::rowLength - 1
    )
    {
        FatalErrorIn
        (
            "direction TensorN<Cmpt, length>::cmpt()\n"
            "(\n"
            "    const direction i,\n"
            "    const direction j\n"
            ") const"
        )   << "Direction out of range (0 "
            << TensorN<Cmpt, length>::rowLength - 1 << ")"
            << " and (i j) = ("  << i << " " << j << ")"
            << abort(FatalError);
    }

    return i*TensorN<Cmpt, length>::rowLength + j;
}


template <class Cmpt, int length>
inline const Cmpt& TensorN<Cmpt, length>::operator()
(
    const direction i,
    const direction j
) const
{
    return this->operator[](i*TensorN<Cmpt, length>::rowLength + j);
}


template <class Cmpt, int length>
inline Cmpt& TensorN<Cmpt, length>::operator()
(
    const direction i,
    const direction j
)
{
    return this->operator[](i*TensorN<Cmpt, length>::rowLength + j);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Inner-product between two tensors
template <class Cmpt, int length>
inline typename
innerProduct<TensorN<Cmpt, length>, TensorN<Cmpt, length> >::type
operator&(const TensorN<Cmpt, length>& t1, const TensorN<Cmpt, length>& t2)
{
    TensorN<Cmpt, length> result(TensorN<Cmpt, length>::zero);

    label i = 0;
    label j = 0;

    label m, n;

    for (label row = 0; row < TensorN<Cmpt, length>::rowLength; row++)
    {
        for (label col = 0; col < TensorN<Cmpt, length>::rowLength; col++)
        {
            Cmpt& r = result.v_[i];
            m = j;
            n = col;

            for
            (
                label row2 = 0;
                row2 < TensorN<Cmpt, length>::rowLength;
                row2++
            )
            {
                r += t1.v_[m]*t2.v_[n];
                m++;
                n += TensorN<Cmpt, length>::rowLength;
            }
            i++;
        }
        j += TensorN<Cmpt, length>::rowLength;
    }

    return result;
}

//- Inner-product between diagonal tensor and tensor
template <class Cmpt, int length>
inline typename
innerProduct<DiagTensorN<Cmpt, length>, TensorN<Cmpt, length> >::type
operator&
(
    const DiagTensorN<Cmpt, length>& dt1,
    const TensorN<Cmpt, length>& t2
)
{
    TensorN<Cmpt, length> result;

    label k = 0;
    for (label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        const Cmpt& xx = dt1.v_[i];

        for (label j = 0; j < TensorN<Cmpt, length>::rowLength; j++)
        {
            result.v_[k] = xx*t2.v_[k];
            k++;
        }
    }

    return result;
}

//- Inner-product between tensor and diagonal tensor
template <class Cmpt, int length>
inline typename
innerProduct<TensorN<Cmpt, length>, DiagTensorN<Cmpt, length> >::type
operator&
(
    const TensorN<Cmpt, length>& t1,
    const DiagTensorN<Cmpt, length>& dt2
)
{
    TensorN<Cmpt, length> result;

    label k = 0;
    for (label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        for (label j = 0; j < TensorN<Cmpt, length>::rowLength; j++)
        {
            result.v_[k] = t1.v_[k]*dt2.v_[j];
            k++;
        }
    }

    return result;
}


//- Inner-product between spherical tensor and tensor
template <class Cmpt, int length>
inline typename
innerProduct<SphericalTensorN<Cmpt, length>, TensorN<Cmpt, length> >::type
operator&
(
    const SphericalTensorN<Cmpt, length>& st1,
    const TensorN<Cmpt, length>& t2
)
{
    const Cmpt& s = st1.v_[0];
    TensorN<Cmpt, length> res;
    VectorSpaceOps<TensorN<Cmpt, length>::nComponents,0>::opSV
    (
        res,
        s,
        t2,
        multiplyOp<Cmpt>()
    );

    return res;
}


//- Inner-product between tensor and spherical tensor
template <class Cmpt, int length>
inline typename
innerProduct<TensorN<Cmpt, length>, SphericalTensorN<Cmpt, length> >::type
operator&
(
    const TensorN<Cmpt, length>& t1,
    const SphericalTensorN<Cmpt, length>& st2
)
{
    const Cmpt& s = st2.v_[0];
    TensorN<Cmpt, length> res;
    VectorSpaceOps<TensorN<Cmpt, length>::nComponents,0>::opVS
    (
        res,
        t1,
        s,
        multiplyOp<Cmpt>()
    );

    return res;
}


//- Inner-product between a tensor and a vector
template <class Cmpt, int length>
inline typename
innerProduct<TensorN<Cmpt, length>, VectorN<Cmpt, length> >::type
operator&(const TensorN<Cmpt, length>& t, const VectorN<Cmpt, length>& v)
{
    VectorN<Cmpt, length> result(VectorN<Cmpt, length>::zero);

    label i = 0;
    for (label row = 0; row < TensorN<Cmpt, length>::rowLength; row++)
    {
        Cmpt& r = result.v_[row];

        for (label col = 0; col < TensorN<Cmpt, length>::rowLength; col++)
        {
            r += t.v_[i]*v.v_[col];
            i++;
        }
    }

    return result;
}


//- Inner-product between a vector and a tensor
template <class Cmpt, int length>
inline typename
innerProduct<VectorN<Cmpt, length>, TensorN<Cmpt, length> >::type
operator&(const VectorN<Cmpt, length>& v, const TensorN<Cmpt, length>& t)
{
    VectorN<Cmpt, length> result(VectorN<Cmpt, length>::zero);

    for (label col = 0; col < TensorN<Cmpt, length>::rowLength; col++)
    {
        label j = col;
        Cmpt& r = result.v_[col];

        for (label row = 0; row < TensorN<Cmpt, length>::rowLength; row++)
        {
            r += v.v_[row]*t.v_[j];
            j += TensorN<Cmpt, length>::rowLength;
        }
    }

    return result;
}


//- Outer-product between two vectors
template <class Cmpt, int length>
inline typename
outerProduct<VectorN<Cmpt, length>, VectorN<Cmpt, length> >::type
operator*(const VectorN<Cmpt, length>& v1, const VectorN<Cmpt, length>& v2)
{
    TensorN<Cmpt, length> result(TensorN<Cmpt, length>::zero);

    label i = 0;
    for (label row = 0; row < TensorN<Cmpt, length>::rowLength; row++)
    {
        for (label col = 0; col < TensorN<Cmpt, length>::rowLength; col++)
        {
            result.v_[i] = v1.v_[row]*v2.v_[col];
            i++;
        }
    }

    return result;
}


//- Addition of TensorN and TensorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length>
operator+(const TensorN<Cmpt, length>& t1, const TensorN<Cmpt, length>& t2)
{
    TensorN<Cmpt, length> res;
    VectorSpaceOps<TensorN<Cmpt, length>::nComponents,0>::op
    (
        res,
        t1,
        t2,
        plusOp<Cmpt>()
    );

    return res;
}


//- Addition of TensorN and DiagTensorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length>
operator+
(
    const TensorN<Cmpt, length>& t1,
    const DiagTensorN<Cmpt, length>& dt2
)
{
    TensorN<Cmpt, length> result(t1);

    label diag = 0;
    for (label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result.v_[diag] += dt2.v_[i];
        diag += TensorN<Cmpt, length>::rowLength + 1;
    }

    return result;
}


//- Addition of DiagTensorN and TensorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length>
operator+
(
    const DiagTensorN<Cmpt, length>& dt1,
    const TensorN<Cmpt, length>& t2
)
{
    TensorN<Cmpt, length> result(t2);

    label diag = 0;
    for (label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result.v_[diag] += dt1.v_[i];
        diag += TensorN<Cmpt, length>::rowLength + 1;
    }

    return result;
}


//- Addition of TensorN and SphericalTensorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length>
operator+
(
    const TensorN<Cmpt, length>& t1,
    const SphericalTensorN<Cmpt, length>& st2
)
{
    TensorN<Cmpt, length> result(t1);

    const Cmpt& s = st2.v_[0];
    label diag = 0;
    for (label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result.v_[diag] += s;
        diag += TensorN<Cmpt, length>::rowLength + 1;
    }

    return result;
}


//- Addition of SphericalTensorN and TensorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length>
operator+
(
    const SphericalTensorN<Cmpt, length>& st1,
    const TensorN<Cmpt, length>& t2
)
{
    TensorN<Cmpt, length> result(t2);

    const Cmpt& s = st1.v_[0];
    label diag = 0;
    for (label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result.v_[diag] += s;
        diag += TensorN<Cmpt, length>::rowLength + 1;
    }

    return result;
}


//- Subtraction of TensorN and TensorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length>
operator-(const TensorN<Cmpt, length>& t1, const TensorN<Cmpt, length>& t2)
{
    TensorN<Cmpt, length> res;
    VectorSpaceOps<TensorN<Cmpt, length>::nComponents,0>::op
    (
        res,
        t1,
        t2,
        minusOp<Cmpt>()
    );

    return res;
}


//- Subtraction of TensorN and DiagTensorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length>
operator-
(
    const TensorN<Cmpt, length>& t1,
    const DiagTensorN<Cmpt, length>& dt2
)
{
    TensorN<Cmpt, length> result(t1);

    label diag = 0;
    for (label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result.v_[diag] -= dt2.v_[i];
        diag += TensorN<Cmpt, length>::rowLength + 1;
    }

    return result;
}


//- Subtraction of DiagTensorN and TensorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length>
operator-
(
    const DiagTensorN<Cmpt, length>& dt1,
    const TensorN<Cmpt, length>& t2
)
{
    TensorN<Cmpt, length> result(-t2);

    label diag = 0;
    for (label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result.v_[diag] += dt1.v_[i];
        diag += TensorN<Cmpt, length>::rowLength + 1;
    }

    return result;
}


//- Subtraction of TensorN and SphericalTensorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length>
operator-
(
    const TensorN<Cmpt, length>& t1,
    const SphericalTensorN<Cmpt, length>& st2
)
{
    TensorN<Cmpt, length> result(t1);

    const Cmpt& s = st2.v_[0];
    label diag = 0;
    for (label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result.v_[diag] -= s;
        diag += TensorN<Cmpt, length>::rowLength + 1;
    }

    return result;
}


//- Subtraction of SphericalTensorN and TensorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length>
operator-
(
    const SphericalTensorN<Cmpt, length>& st1,
    const TensorN<Cmpt, length>& t2
)
{
    TensorN<Cmpt, length> result(-t2);

    const Cmpt& s = st1.v_[0];
    label diag = 0;
    for (label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result.v_[diag] += s;
        diag += TensorN<Cmpt, length>::rowLength + 1;
    }

    return result;
}


//- Division of a scalar by a TensorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length>
operator/(const scalar s, const TensorN<Cmpt, length>& t)
{
    return s*inv(t);
}

//- Inner Product of a VectorN by an inverse TensorN
template <class Cmpt, int length>
inline VectorN<Cmpt, length>
operator/(const VectorN<Cmpt, length>& v, const TensorN<Cmpt, length>& t)
{
    return v & inv(t);
}

//- Inner Product of a TensorN by an inverse TensorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length>
operator/(const TensorN<Cmpt, length>& t1, const TensorN<Cmpt, length>& t2)
{
    return t1 & inv(t2);
}


//- Inner Product of a DiagTensorN and an inverse TensorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length>
operator/
(
    const DiagTensorN<Cmpt, length>& dt1,
    const TensorN<Cmpt, length>& t2
)
{
    return dt1 & inv(t2);
}


//- Inner Product of a TensorN and an inverse DiagTensorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length>
operator/
(
    const TensorN<Cmpt, length>& t1,
    const DiagTensorN<Cmpt, length>& dt2
)
{
    return t1 & inv(dt2);
}


//- Inner Product of a SphericalTensorN and an inverse TensorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length>
operator/
(
    const SphericalTensorN<Cmpt, length>& st1,
    const TensorN<Cmpt, length>& t2
)
{
    return st1.v_[0] * inv(t2);
}


//- Inner Product of a TensorN and an inverse SphericalTensorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length>
operator/
(
    const TensorN<Cmpt, length>& t1,
    const SphericalTensorN<Cmpt, length>& st2
)
{
    TensorN<Cmpt, length> result;

    const Cmpt& s = st2[0];
    for (label i = 0; i < TensorN<Cmpt, length>::nComponents; i++)
    {
        result.v_[i] = t1.v_[i]/s;
    }

    return result;
}


template <class Cmpt, int length>
inline Cmpt det(const TensorN<Cmpt, length>& t)
{
    // Calculate determinant via sub-determinants
    Cmpt result = pTraits<Cmpt>::zero;

    TensorN<Cmpt, length - 1> subMatrix;

    // Laplace expansion using first column
    // Bugfix, TU, 31 May 2016

    for (label i = 0; i < length; i++)
    {
        // Column of submatrix
        label js = 0;

        for (label j = 1; j < length; j++)
        {
            // Row of submatrix
            label is = 0;

            for (label k = 0; k < length; k++)
            {
                // Row which is crossed out in the expansion
                if (k == i)
                {
                    continue;
                }
                else
                {
                    subMatrix(is,js) = t(k,j);
                    is++;
                }
            }
            js++;
        }

        // Handle +/- sign switch
        result += pow(-1, i)*t(i, 0)*det(subMatrix);
    }

    return result;
}


// Determinant: partial specialisation for rank 1
template<class Cmpt>
inline Cmpt det(const TensorN<Cmpt, 1>& t)
{
    return t(0, 0);
}


// Determinant: partial specialisation for rank 2
template<class Cmpt>
inline Cmpt det(const TensorN<Cmpt, 2>& t)
{
    return t(0, 0)*t(1, 1) - t(0, 1)*t(1, 0);
}


//- Return the inverse of a tensor give the determinant
//  Uses Gauss-Jordan Elimination with full pivoting
template <class Cmpt, int length>
inline TensorN<Cmpt, length> inv(const TensorN<Cmpt, length>& t)
{
    TensorN<Cmpt, length> result(t);

    label iRow = 0, iCol = 0;
    Cmpt largestCoeff, temp;
    Cmpt* __restrict__ srcIter;
    Cmpt* __restrict__ destIter;

    // Lists used for bookkeeping on the pivoting
    List<label> indexCol(length), indexRow(length), iPivot(length);

    iPivot = 0;

    label curRowOffset;

    // Main loop over columns to be reduced
    for (label i = 0; i < length; i++)
    {
        largestCoeff = pTraits<Cmpt>::zero;

        // Search for pivot element
        curRowOffset = 0;

        for (label j = 0; j < length; j++)
        {
            if (iPivot[j] != 1)
            {
                for (int k = 0; k < length; k++)
                {
                    if (iPivot[k] == 0)
                    {
                        if
                        (
                            (temp = Foam::mag(result[curRowOffset+k]))
                         >= largestCoeff
                        )
                        {
                            largestCoeff = temp;
                            iRow = j;
                            iCol = k;
                        }
                    }
                }
            }
            curRowOffset += length;
        }
        ++(iPivot[iCol]);

        // We now have the pivot element
        // Interchange rows if needed
        if (iRow != iCol)
        {
            srcIter = &result(iRow, 0);
            destIter = &result(iCol, 0);

            for (label j = 0; j < length; j++)
            {
                Swap((*srcIter), (*destIter));
                srcIter++;
                destIter++;
            }
        }
        indexRow[i] = iRow;
        indexCol[i] = iCol;

        //Check for singularity
        srcIter = &result(iCol, iCol);  // Dummy pointer to reduce indexing
        if ((*srcIter) == Cmpt(0))
        {
            FatalErrorIn
            (
                "inline TensorN<Cmpt, length> inv("
                "const TensorN<Cmpt, length>& t)"
            )   << "Singular tensor" << length << Foam::abort(FatalError);
        }

        // Divide the pivot row by pivot element
        temp = pTraits<Cmpt>::one/(*srcIter);
        (*srcIter) = pTraits<Cmpt>::one;

        srcIter = &result(iCol, 0);
        for (label j = 0; j < length; j++)
        {
            (*srcIter) *= temp;
            srcIter++;
        }

        // Reduce the rows, excluding the pivot row
        for (label j = 0; j < length; j++)
        {
            if (j != iCol)
            {
                destIter = &result(j, 0);
                srcIter = &result(iCol, 0);

                temp=destIter[iCol];
                destIter[iCol] = pTraits<Cmpt>::zero;

                for (label k = 0; k < length; k++)
                {
                    (*destIter) -= (*srcIter)*temp;
                    srcIter++;
                    destIter++;
                }
            }
        }
    }

    // Unscamble the solution
    for (label i = length - 1; i >= 0; i--)
    {
        if (indexRow[i] != indexCol[i])
        {
            srcIter = &result[indexRow[i]];
            destIter = &result[indexCol[i]];
            for (label j = 0; j < length; j++)
            {
                Swap((*srcIter), (*destIter));
                srcIter += length;
                destIter += length;
            }
        }
    }

    return result;
}


//- Inverse: partial template specialisation for rank 1
template <class Cmpt>
inline TensorN<Cmpt, 1> inv(const TensorN<Cmpt, 1>& t)
{
    return TensorN<Cmpt, 1>(1/t(0, 0));
}


//- Inverse: partial template specialisation for rank 2
template <class Cmpt>
inline TensorN<Cmpt, 2> inv(const TensorN<Cmpt, 2>& t)
{
    TensorN<Cmpt, 2> result(t);

    Cmpt oneOverDet = 1/(t(0, 0)*t(1, 1) - t(0, 1)*t(1, 0));

    result(0, 0) = oneOverDet*t(1, 1);
    result(0, 1) = -oneOverDet*t(0, 1);
    result(1, 0) = -oneOverDet*t(1, 0);
    result(1, 1) = oneOverDet*t(0, 0);

    return result;
}


//- Return tensor diagonal
template <class Cmpt, int length>
inline DiagTensorN<Cmpt, length> diag(const TensorN<Cmpt, length>& t)
{
    return t.diag();
}


//- Return tensor diagonal
template <class Cmpt, int length>
inline TensorN<Cmpt, length> negSumDiag(const TensorN<Cmpt, length>& t)
{
    return t.negSumDiag();
}


//- Scale tensor row with a VectorN
template <class Cmpt, int length>
inline TensorN<Cmpt, length> scaleRow
(
    const TensorN<Cmpt, length>& t,
    const VectorN<Cmpt, length>& v
)
{
    TensorN<Cmpt, length> result;

    for (label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        for (label j = 0; j < TensorN<Cmpt, length>::rowLength; j++)
        {
            result(i, j) = t(i, j)*v(j);
        }
    }

    return result;
}


//- Return the component sum
// template <class Cmpt, int length>
// inline Cmpt sum(const TensorN<Cmpt, length>& t)
// {
//     Cmpt result=Cmpt::zero;
//     for(int i=0; i<TensorN<Cmpt, length>::nComponents; i++)
//     {
//         result += t[i];
//     }
//     return result;
// }

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Cmpt, int length>
class outerProduct<TensorN<Cmpt, length>, Cmpt>
{
public:

    typedef TensorN<Cmpt, length> type;
};

template<class Cmpt, int length>
class outerProduct<Cmpt, TensorN<Cmpt, length> >
{
public:

    typedef TensorN<Cmpt, length> type;
};

template<class Cmpt, int length>
class innerProduct<DiagTensorN<Cmpt, length>, TensorN<Cmpt, length> >
{
public:

    typedef TensorN<Cmpt, length> type;
};


template<class Cmpt, int length>
class innerProduct<TensorN<Cmpt, length>, DiagTensorN<Cmpt, length> >
{
public:

    typedef TensorN<Cmpt, length> type;
};


template<class Cmpt, int length>
class innerProduct<SphericalTensorN<Cmpt, length>, TensorN<Cmpt, length> >
{
public:

    typedef TensorN<Cmpt, length> type;
};


template<class Cmpt, int length>
class innerProduct<TensorN<Cmpt, length>, SphericalTensorN<Cmpt, length> >
{
public:

    typedef TensorN<Cmpt, length> type;
};


template<class Cmpt, int length>
class innerProduct<VectorN<Cmpt, length>, TensorN<Cmpt, length> >
{
public:

    typedef VectorN<Cmpt, length> type;
};


template<class Cmpt, int length>
class innerProduct<TensorN<Cmpt, length>, VectorN<Cmpt, length> >
{
public:

    typedef VectorN<Cmpt, length> type;
};


template<class Cmpt, int length>
class innerProduct<TensorN<Cmpt, length>, TensorN<Cmpt, length> >
{
public:

    typedef TensorN<Cmpt, length> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
