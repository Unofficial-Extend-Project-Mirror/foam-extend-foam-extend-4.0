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

#include "BlockCoeff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

template<class Type>
tmp<CoeffField<Type> > inv(const CoeffField<Type>& f)
{
    // The inverse of a linear coefficient type is currently done "by
    // hand".  The need for this will disappear once the diagonal tensor
    // type is introduced.  HJ, 24/May/2005

    typedef typename CoeffField<Type>::linearTypeField fieldType;
    typedef typename CoeffField<Type>::linearType linearType;

    // Create result
    tmp<CoeffField<Type> > tresult(new CoeffField<Type>(f.size()));
    CoeffField<Type>& result = tresult();

    if (f.activeType() == blockCoeffBase::SCALAR)
    {
        result = 1.0/f.asScalar();
    }
    else if (f.activeType() == blockCoeffBase::LINEAR)
    {
        const fieldType& lf = f.asLinear();

        fieldType inverse =
            cmptDivide
            (
                fieldType(lf.size(), pTraits<linearType>::one),
                lf
            );

        result = inverse;
    }
    else if (f.activeType() == blockCoeffBase::SQUARE)
    {
        result = inv(f.asSquare());
    }

    return tresult;
}


template<class Type>
void multiply
(
    Field<Type>& f,
    const CoeffField<Type>& f1,
    const Type& f2
)
{
    if (f1.activeType() == blockCoeffBase::SCALAR)
    {
        f = f1.asScalar()*f2;
    }
    else if (f1.activeType() == blockCoeffBase::LINEAR)
    {
        f = cmptMultiply(f1.asLinear(), f2);
    }
    else if (f1.activeType() == blockCoeffBase::SQUARE)
    {
        f = f1.asSquare() & f2;
    }
}


template<class Type>
void multiply
(
    Field<Type>& f,
    const CoeffField<Type>& f1,
    const Field<Type>& f2
)
{
    if (f1.activeType() == blockCoeffBase::SCALAR)
    {
        f = f1.asScalar()*f2;
    }
    else if (f1.activeType() == blockCoeffBase::LINEAR)
    {
        f = cmptMultiply(f1.asLinear(), f2);
    }
    else if (f1.activeType() == blockCoeffBase::SQUARE)
    {
        f = f1.asSquare() & f2;
    }
}


template<class Type>
void multiply
(
    Field<Type>& f,
    const Field<Type>& f1,
    const CoeffField<Type>& f2
)
{
    if (f2.activeType() == blockCoeffBase::SCALAR)
    {
        f = f1*f2.asScalar();
    }
    else if (f2.activeType() == blockCoeffBase::LINEAR)
    {
        f = cmptMultiply(f1, f2.asLinear());
    }
    else if (f2.activeType() == blockCoeffBase::SQUARE)
    {
        f = f1 & f2.asSquare();
    }
}


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

#define UNARY_OPERATOR(op, opFunc)                                            \
                                                                              \
template<class Type>                                                          \
void opFunc                                                                   \
(                                                                             \
    CoeffField<Type>& f,                                                      \
    const CoeffField<Type>& f1                                                \
)                                                                             \
{                                                                             \
    typedef CoeffField<Type> TypeCoeffField;                                  \
                                                                              \
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;         \
    typedef typename TypeCoeffField::linearTypeField linearTypeField;         \
    typedef typename TypeCoeffField::squareTypeField squareTypeField;         \
                                                                              \
    if (f.activeType() == blockCoeffBase::SCALAR)                             \
    {                                                                         \
        scalarTypeField sf = f1.asScalar();                                   \
        sf.opFunc();                                                          \
        f = sf;                                                               \
    }                                                                         \
    else if (f.activeType() == blockCoeffBase::LINEAR)                        \
    {                                                                         \
        linearTypeField sf = f1.asLinear();                                   \
        sf.opFunc();                                                          \
        f = sf;                                                               \
    }                                                                         \
    else if (f.activeType() == blockCoeffBase::SQUARE)                        \
    {                                                                         \
        squareTypeField sf = f1.asSquare();                                   \
        sf.opFunc();                                                          \
        f = sf;                                                               \
    }                                                                         \
}                                                                             \
                                                                              \
template<class Type>                                                          \
tmp<CoeffField<Type> > operator op                                            \
(                                                                             \
    const CoeffField<Type>& f1                                                \
)                                                                             \
{                                                                             \
    tmp<CoeffField<Type> > tf(new CoeffField<Type>(f1.size()));               \
    opFunc(tf(), f1);                                                         \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<class Type>                                                          \
tmp<CoeffField<Type> > operator op                                            \
(                                                                             \
    const tmp<CoeffField<Type> >& tf1                                         \
)                                                                             \
{                                                                             \
    tmp<CoeffField<Type> > tf(tf1.ptr());                                     \
    opFunc(tf(), tf());                                                       \
    return tf;                                                                \
}

UNARY_OPERATOR(-, negate)

#undef UNARY_OPERATOR


#define BINARY_OPERATOR_FF(Type1, Type2, op, opFunc)                          \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const CoeffField<Type1>& f1,                                              \
    const Type2& f2                                                           \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(new Field<Type>(f1.size()));                         \
    opFunc(tf(), f1, f2);                                                     \
    return tf;                                                                \
}                                                                             \
                                                                              \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const CoeffField<Type1>& f1,                                              \
    const Field<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(new Field<Type>(f1.size()));                         \
    opFunc(tf(), f1, f2);                                                     \
    return tf;                                                                \
}                                                                             \
                                                                              \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const Field<Type2>& f1,                                                   \
    const CoeffField<Type1>& f2                                               \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(new Field<Type>(f1.size()));                         \
    opFunc(tf(), f1, f2);                                                     \
    return tf;                                                                \
}

#define BINARY_OPERATOR_FTR(Type1, Type2, op, opFunc)                         \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const CoeffField<Type1>& f1,                                              \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(tf2.ptr());                                          \
    opFunc(tf(), f1, tf());                                                   \
    return tf;                                                                \
}

#define BINARY_OPERATOR_FT(Type1, Type2, op, opFunc)                          \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const Field<Type1>& f1,                                                   \
    const tmp<CoeffField<Type2> >& tf2                                        \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf = f1 op tf2();                                       \
    tf2.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_OPERATOR_TRF(Type1, Type2, op, opFunc)                         \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const tmp<CoeffField<Type1> >& tf1,                                       \
    const Field<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(tf1.ptr());                                          \
    opFunc(tf(), tf(), f2);                                                   \
    return tf;                                                                \
}

#define BINARY_OPERATOR_TF(Type1, Type2, op, opFunc)                          \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const tmp<CoeffField<Type1> >& tf1,                                       \
    const Field<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf = tf1() op f2;                                       \
    tf1.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_OPERATOR_TRT(Type1, Type2, op, opFunc)                         \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const tmp<CoeffField<Type1> >& tf1,                                       \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(tf1.ptr());                                          \
    opFunc(tf(), tf(), tf2());                                                \
    tf2.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_OPERATOR_TTR(Type1, Type2, op, opFunc)                         \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<CoeffField<Type2> >& tf2                                        \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(tf2.ptr());                                          \
    opFunc(tf(), tf1(), tf());                                                \
    tf1.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_OPERATOR_R(Type1, Type2, op, opFunc)                           \
    BINARY_OPERATOR_FF(Type1, Type2, op, opFunc)                              \
    BINARY_OPERATOR_FTR(Type1, Type2, op, opFunc)                             \
    BINARY_OPERATOR_TRF(Type1, Type2, op, opFunc)                             \
    BINARY_OPERATOR_TRT(Type1, Type2, op, opFunc)

// Operator multiply is not available for all types, as it expands rank
// HJ, 17/Jun/2010
BINARY_OPERATOR_R(Type, Type, *, multiply)


#undef BINARY_OPERATOR_R
#undef BINARY_OPERATOR_FF
#undef BINARY_OPERATOR_FTR
#undef BINARY_OPERATOR_TF
#undef BINARY_OPERATOR_TTR
#undef BINARY_OPERATOR_FT
#undef BINARY_OPERATOR_TRF
#undef BINARY_OPERATOR_TRT


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
