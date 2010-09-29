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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template <class Cmpt, int length>
const char* const TensorN<Cmpt, length>::typeName =
    string("tensor" + name(length)).c_str();

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
    for (register label i = 0; i < TensorN<Cmpt, length>::nComponents; i++)
    {
        this->v_[i] = tx;
    }
}


//- Construct from Istream
template <class Cmpt, int length>
inline TensorN<Cmpt, length>::TensorN(Istream& is)
:
    VectorSpace<TensorN<Cmpt, length>, Cmpt, length*length>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return tensor transpose
template <class Cmpt, int length>
inline TensorN<Cmpt, length> TensorN<Cmpt, length>::T() const
{
    TensorN<Cmpt, length> transpose;

    for (register label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        for (register label j = 0; j < TensorN<Cmpt, length>::rowLength; j++)
        {
            transpose(i, j) = this->operator()(j, i);
        }
    }

    return transpose;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

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

    for (register label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        for (register label j = 0; j < TensorN<Cmpt, length>::rowLength; j++)
        {
            Cmpt& r = result(i, j);

            for (int k = 0; k < TensorN<Cmpt, length>::rowLength; k++) 
            {
                r += t1(i, k)*t2(k, j);
            }
        }
    }

    return result;
}


//- Inner-product between a tensor and a vector
template <class Cmpt, int length>
inline typename
innerProduct<TensorN<Cmpt, length>, VectorN<Cmpt, length> >::type
operator&(const TensorN<Cmpt, length>& t, const VectorN<Cmpt, length>& v)
{
    VectorN<Cmpt, length> result(VectorN<Cmpt, length>::zero);

    for (register label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        Cmpt& r = result(i);

        for (register label j = 0; j < TensorN<Cmpt, length>::rowLength; j++)
        {
            r += t(i, j)*v(j);
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

    for (register label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        Cmpt& r = result(i);

        for (register label j = 0; j < TensorN<Cmpt, length>::rowLength; j++)
        {
            r += v(i)*t(j, i);
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

    for (register label i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        for (register label j = 0; j < TensorN<Cmpt, length>::rowLength; j++)
        {
            result(i, j) = v1(i)*v2(j);
        }
    }

    return result;
}


// UNOPTIMIZED VERSION
/*
//- Return the inverse of a tensor give the determinant
//  Uses Gauss-Jordan Elimination with full pivoting
template <class Cmpt, int length>
inline TensorN<Cmpt, length> inv(const TensorN<Cmpt, length>& t)
{
    TensorN<Cmpt, length> result(t);

    label i, j, k, iRow=0, iCol=0;
    Cmpt bigValue, temp, pivotInv;

    // Lists used for bookkeeping on the pivoting
    List<label> indexCol(length), indexRow(length), iPivot(length);

    iPivot=0;

    // Main loop over columns to be reduced
    for (i=0; i<length; i++)
    {
        bigValue = pTraits<Cmpt>::zero;

        //Search for pivot element
        for (j=0; j<length; j++)
        {
            if (iPivot[j] != 1)
            {
                for (k=0; k<length; k++)
                {
                    if (iPivot[k] == 0)
                    {
                        if (Foam::mag(result(j,k)) >= bigValue)
                        {
                            bigValue = Foam::mag(result(j,k));
                            iRow = j;
                            iCol = k;
                        }
                    }
                }
            }
        }
        ++(iPivot[iCol]);

        // We now have the pivot element
        // Interchange rows if needed
        if (iRow != iCol)
        {
            for (j=0; j<length; j++)
            {
                Swap(result(iRow,j), result(iCol,j));
            }
        }
        indexRow[i] = iRow;
        indexCol[i] = iCol;

        //Check for singularity
        if (result(iCol, iCol) == 0.0)
        {
            FatalErrorIn("inline TensorN<Cmpt, length> inv(const TensorN<Cmpt, length>& t)")
                << "Singular tensor" << length << Foam::abort(FatalError);
        }

        // Divide the pivot row by pivot element
        pivotInv = pTraits<Cmpt>::one/result(iCol, iCol);
        result(iCol, iCol) = pTraits<Cmpt>::one;

        // Multiply all row elements by inverse
        for (j=0; j<length; j++)
        {
            result(iCol,j) *= pivotInv;
        }

        // Reduce the rows
        for (j=0; j<length; j++)
        {
            if (j != iCol)
            {
                temp=result(j,iCol);
                result(j,iCol) = pTraits<Cmpt>::zero;

                for (k=0; k<length; k++)
                {
                    result(j,k) -= result(iCol,k)*temp;
                }
            }
        }
    }

    // Unscamble the solution
    for (i=length-1; i>=0; i--)
    {
        if (indexRow[i] != indexCol[i])
        {
            for (j=0; j<length; j++)
            {
                Swap(result(j,indexRow[i]), result(j,indexCol[i]));
            }
        }
    }

    return result;
}
*/

//- Return the inverse of a tensor give the determinant
//  Uses Gauss-Jordan Elimination with full pivoting
template <class Cmpt, int length>
inline TensorN<Cmpt, length> inv(const TensorN<Cmpt, length>& t)
{
    TensorN<Cmpt, length> result(t);

    label iRow=0, iCol=0;
    Cmpt largestCoeff, temp;
    Cmpt *srcIter, *destIter;

    // Lists used for bookkeeping on the pivoting
    List<label> indexCol(length), indexRow(length), iPivot(length);

    iPivot=0;

    // Main loop over columns to be reduced
    for (register label i=0; i<length; i++)
    {
        largestCoeff = pTraits<Cmpt>::zero;

        //Search for pivot element
        register label curRowOffset = 0;
        for (register label j = 0; j < length; j++)
        {
            if (iPivot[j] != 1)
            {
                for (register label k = 0; k < length; k++)
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
            srcIter = &result(iRow,0);
            destIter = &result(iCol,0);

            for (register label j=0; j<length; j++)
            {
                Swap((*srcIter), (*destIter));
                srcIter++;
                destIter++;
            }
        }
        indexRow[i] = iRow;
        indexCol[i] = iCol;

        //Check for singularity
        srcIter = &result(iCol, iCol);  //Dummy pointer to reduce indexing
        if ((*srcIter) == Cmpt(0.0))
        {
            FatalErrorIn("inline TensorN<Cmpt, length> inv(const TensorN<Cmpt, length>& t)")
                << "Singular tensor" << length << Foam::abort(FatalError);
        }

        // Divide the pivot row by pivot element
        temp = pTraits<Cmpt>::one/(*srcIter);
        (*srcIter) = pTraits<Cmpt>::one;

        srcIter = &result(iCol,0);
        for (register label j=0; j<length; j++)
        {
            (*srcIter) *= temp;
            srcIter++;
        }

        // Reduce the rows, excluding the pivot row
        for (register label j=0; j<length; j++)
        {
            if (j != iCol)
            {
                destIter = &result(j,0);
                srcIter = &result(iCol,0);

                temp=destIter[iCol];
                destIter[iCol] = pTraits<Cmpt>::zero;

                for (register label k=0; k<length; k++)
                {
                    (*destIter) -= (*srcIter)*temp;
                    srcIter++;
                    destIter++;
                }
            }
        }
    }

    // Unscamble the solution
    for (register label i=length-1; i>=0; i--)
    {
        if (indexRow[i] != indexCol[i])
        {
            srcIter = &result[indexRow[i]];
            destIter = &result[indexCol[i]];
            for (register label j=0; j<length; j++)
            {
                Swap((*srcIter), (*destIter));
                srcIter += length;
                destIter += length;
            }
        }
    }

    return result;
}

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

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
