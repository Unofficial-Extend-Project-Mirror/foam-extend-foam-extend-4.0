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

Description
    Generic coefficient field type.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

template<class Type>
tmp<DecoupledCoeffField<Type> > inv(const DecoupledCoeffField<Type>& f);


template<class Type>
void multiply
(
    Field<Type>& f,
    const DecoupledCoeffField<Type>& f1,
    const Type& f2
);


template<class Type>
void multiply
(
    Field<Type>& f,
    const DecoupledCoeffField<Type>& f1,
    const Field<Type>& f2
);


template<class Type>
void multiply
(
    Field<Type>& f,
    const Field<Type>& f1,
    const DecoupledCoeffField<Type>& f2
);


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

#define UNARY_OPERATOR(op, opFunc)                                            \
                                                                              \
template<class Type>                                                          \
void opFunc                                                                   \
(                                                                             \
    DecoupledCoeffField<Type>& f,                                             \
    const DecoupledCoeffField<Type>& f1                                       \
);                                                                            \
                                                                              \
template<class Type>                                                          \
tmp<DecoupledCoeffField<Type> > operator op                                   \
(                                                                             \
    const DecoupledCoeffField<Type>& f1                                       \
);                                                                            \
                                                                              \
template<class Type>                                                          \
tmp<DecoupledCoeffField<Type> > operator op                                   \
(                                                                             \
    const tmp<DecoupledCoeffField<Type> >& tf1                                \
);

UNARY_OPERATOR(-, negate)

#undef UNARY_OPERATOR


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

#define BINARY_OPERATOR_FF(Type1, Type2, op, opFunc)                          \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const DecoupledCoeffField<Type1>& f1,                                     \
    const Type2& f2                                                           \
);                                                                            \
                                                                              \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const DecoupledCoeffField<Type1>& f1,                                     \
    const Field<Type2>& f2                                                    \
);                                                                            \
                                                                              \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const Field<Type2>& f1,                                                   \
    const DecoupledCoeffField<Type1>& f2                                      \
);


#define BINARY_OPERATOR_FTR(Type1, Type2, op, opFunc)                         \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const DecoupledCoeffField<Type1>& f1,                                     \
    const tmp<Field<Type2> >& tf2                                             \
);

#define BINARY_OPERATOR_FT(Type1, Type2, op, opFunc)                          \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const Field<Type1>& f1,                                                   \
    const tmp<DecoupledCoeffField<Type2> >& tf2                               \
);

#define BINARY_OPERATOR_TRF(Type1, Type2, op, opFunc)                         \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const DecoupledCoeffField<Type2>& f2                                      \
);

#define BINARY_OPERATOR_TF(Type1, Type2, op, opFunc)                          \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const tmp<DecoupledCoeffField<Type1> >& tf1,                              \
    const Field<Type2>& f2                                                    \
);

#define BINARY_OPERATOR_TRT(Type1, Type2, op, opFunc)                         \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const tmp<DecoupledCoeffField<Type1> >& tf1,                              \
    const tmp<Field<Type2> >& tf2                                             \
);

#define BINARY_OPERATOR_TTR(Type1, Type2, op, opFunc)                         \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<DecoupledCoeffFieldy<Type2> >& tf2                              \
);

#define BINARY_OPERATOR_R(Type1, Type2, op, opFunc)                           \
    BINARY_OPERATOR_FF(Type1, Type2, op, opFunc)                              \
    BINARY_OPERATOR_FTR(Type1, Type2, op, opFunc)                             \
    BINARY_OPERATOR_TRF(Type1, Type2, op, opFunc)                             \
    BINARY_OPERATOR_TRT(Type1, Type2, op, opFunc)

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
