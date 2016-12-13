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

Class
    ExpandTensorN

Description

Author

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return the average of a vector of as a scalar
template <class Cmpt, int length>
inline void contractScalar(Cmpt& result, const VectorN<Cmpt, length>& t)
{
    result = pTraits<Cmpt>::zero;

    for (register int i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result += t[i];
    }

    // Modified 2009/11/3 by I. Clifford
    result *= 1.0/VectorN<Cmpt, length>::nComponents;
}


//- Return the Average of a vector of as a scalar
template <class Cmpt, int length>
inline Cmpt contractScalar(const VectorN<Cmpt, length>& t)
{
    Cmpt result;
    contractScalar(result, t);
    return result;
}


//- Return the diagonal of a TensorN as a scalar
template <class Cmpt, int length>
inline void contractScalar(Cmpt& result, const TensorN<Cmpt, length>& t)
{
    result = pTraits<Cmpt>::zero;

    register int j = 0;
    for (register int i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result += t[j];
        j += TensorN<Cmpt, length>::rowLength + 1;
    }

    // Modified 2009/11/3 by I. Clifford
    result *= 1.0/TensorN<Cmpt, length>::rowLength;
}


//- Return the diagonal of a TensorN as a scalar
template <class Cmpt, int length>
inline Cmpt contractScalar(const TensorN<Cmpt, length>& t)
{
    Cmpt result;
    contractScalar(result, t);
    return result;
}


//- Return the diagonal of a DiagTensorN as a scalar
template <class Cmpt, int length>
inline void contractScalar(Cmpt& result, const DiagTensorN<Cmpt, length>& t)
{
    result = pTraits<Cmpt>::zero;

    for (register int i = 0; i < DiagTensorN<Cmpt, length>::rowLength; i++)
    {
        result += t[i];
    }

    result *= 1.0/TensorN<Cmpt, length>::rowLength;
}


//- Return the diagonal of a DiagTensorN as a scalar
template <class Cmpt, int length>
inline Cmpt contractScalar(const DiagTensorN<Cmpt, length>& t)
{
    Cmpt result;
    contractScalar(result, t);
    return result;
}


//- Return the diagonal of a SphericalTensorN as a scalar
template <class Cmpt, int length>
inline void contractScalar
(
    Cmpt& result,
    const SphericalTensorN<Cmpt, length>& t
)
{
    result = t[0];
}


//- Return the diagonal of a SphericalTensorN as a scalar
template <class Cmpt, int length>
inline Cmpt contractScalar(const SphericalTensorN<Cmpt, length>& t)
{
    Cmpt result;
    contractScalar(result, t);
    return result;
}


//- Return the diagonal of a TensorN as a vectorN
template <class Cmpt, int length>
inline void contractLinear
(
    VectorN<Cmpt, length>& result,
    const TensorN<Cmpt, length>& t
)
{
    register int j = 0;
    for (register int i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result[i] = t[j];
        j += TensorN<Cmpt, length>::rowLength + 1;
    }
}


//- Return the diagonal of a TensorN as a vectorN
template <class Cmpt, int length>
inline VectorN<Cmpt, length> contractLinear(const TensorN<Cmpt, length>& t)
{
    VectorN<Cmpt, length> result;
    contractLinear(result, t);
    return result;
}


//- Return the diagonal of a DiagTensorN as a vectorN
template <class Cmpt, int length>
inline void contractLinear
(
    VectorN<Cmpt, length>& result,
    const DiagTensorN<Cmpt, length>& t
)
{
    for (register int i = 0; i < DiagTensorN<Cmpt, length>::rowLength; i++)
    {
        result[i] = t[i];
    }
}


//- Return the diagonal of a DiagTensorN as a vectorN
template <class Cmpt, int length>
inline VectorN<Cmpt, length> contractLinear(const DiagTensorN<Cmpt, length>& t)
{
    VectorN<Cmpt, length> result;
    contractLinear(result, t);
    return result;
}


//- Return the diagonal of a SphericalTensorN as a vectorN
template <class Cmpt, int length>
inline void contractLinear
(
    VectorN<Cmpt, length>& result,
    const SphericalTensorN<Cmpt, length>& t
)
{
    for
    (
        register int i = 0;
        i < SphericalTensorN<Cmpt, length>::rowLength;
        i++
    )
    {
        result[i] = t[0];
    }
}


//- Return the diagonal of a SphericalTensorN as a vectorN
template <class Cmpt, int length>
inline VectorN<Cmpt, length> contractLinear
(
    const SphericalTensorN<Cmpt, length>& t
)
{
    VectorN<Cmpt, length> result;
    contractLinear(result, t);
    return result;
}


//- Return the VectorN given a scalar
template <class Cmpt, int length>
inline void expandScalar(VectorN<Cmpt, length>& result, const Cmpt& v)
{
    for (register int i = 0; i < VectorN<Cmpt, length>::nComponents; i++)
    {
        result[i] = v;
    }
}


//- Return the TensorN given a scalar
template <class Cmpt, int length>
inline void expandScalar(TensorN<Cmpt, length>& result, const Cmpt& v)
{
    result = TensorN<Cmpt, length>::zero;

    register int j = 0;
    for (register int i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result[j] = v;
        j += TensorN<Cmpt, length>::rowLength + 1;
    }
}


//- Return the DiagTensorN given a scalar
template <class Cmpt, int length>
inline void expandScalar(DiagTensorN<Cmpt, length>& result, const Cmpt& v)
{
    for (register int i = 0; i < DiagTensorN<Cmpt, length>::rowLength; i++)
    {
        result[i] = v;
    }
}


//- Return the SphericalTensorN given a scalar
template <class Cmpt, int length>
inline void expandScalar(SphericalTensorN<Cmpt, length>& result, const Cmpt& v)
{
    result[0] = v;
}


//- Return the TensorN given a diagonal vectorN
template <class Cmpt, int length>
inline void expandLinear
(
    TensorN<Cmpt, length>& result,
    const VectorN<Cmpt, length>& v
)
{
    result = TensorN<Cmpt, length>::zero;

    register int j = 0;
    for (register int i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result[j] = v[i];
        j += TensorN<Cmpt, length>::rowLength + 1;
    }
}

//- Return the DiagTensorN given a diagonal vectorN
template <class Cmpt, int length>
inline void expandLinear
(
    DiagTensorN<Cmpt, length>& result,
    const VectorN<Cmpt, length>& v
)
{
    for (register int i = 0; i < DiagTensorN<Cmpt, length>::rowLength; i++)
    {
        result[i] = v[i];
    }
}


//- Return the SphericalTensorN given a diagonal vectorN
template <class Cmpt, int length>
inline void expandLinear
(
    SphericalTensorN<Cmpt, length>& result,
    const VectorN<Cmpt, length>& v
)
{
    result[0] = pTraits<Cmpt>::zero;

    for (register int i = 0; i < DiagTensorN<Cmpt, length>::rowLength; i++)
    {
        result[0] += v[i];
    }

    result[0] *= 1.0/DiagTensorN<Cmpt, length>::rowLength;
}


//- Sum row elements of a TensorN as a vectorN
template <class Cmpt, int length>
inline void sumToDiag
(
    VectorN<Cmpt, length>& result,
    const TensorN<Cmpt, length>& t
)
{
    result = VectorN<Cmpt, length>::zero;

    for (register int i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        for (register int j = 0; j < TensorN<Cmpt, length>::rowLength; j++)
        {
            result[i] += t(i, j);
        }
    }
}


//- Sum row elements of a TensorN as a vectorN
template <class Cmpt, int length>
inline VectorN<Cmpt, length> sumToDiag(const TensorN<Cmpt, length>& t)
{
    VectorN<Cmpt, length> result;
    sumToDiag(result, t);
    return result;
}


//- Sum magnitudes of row elements of a TensorN as a vectorN
template <class Cmpt, int length>
inline void sumMagToDiag
(
    VectorN<Cmpt, length>& result,
    const TensorN<Cmpt, length>& t
)
{
    result = VectorN<Cmpt, length>::zero;

    for (register int i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        for (register int j = 0; j < TensorN<Cmpt, length>::rowLength; j++)
        {
            result[i] += Foam::mag(t(i, j));
        }
    }
}


//- Sum magnitudes of row elements of a TensorN as a vectorN
template <class Cmpt, int length>
inline VectorN<Cmpt, length> sumMagToDiag(const TensorN<Cmpt, length>& t)
{
    VectorN<Cmpt, length> result;
    sumMagToDiag(result, t);
    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
