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
    High performance macro functions for Field\<Type\> algebra.
    These expand using either array element access (for vector machines)
    or pointer dereferencing for scalar machines as appropriate.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(ReturnType, Type1, Func)                               \
                                                                              \
TEMPLATE                                                                      \
void Func                                                                     \
(                                                                             \
    FieldField<Field, ReturnType>& res,                                       \
    const FieldField<Field, Type1>& f                                         \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > Func                                      \
(                                                                             \
    const FieldField<Field, Type1>& f                                         \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > Func                                      \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf                                  \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc)                         \
                                                                              \
TEMPLATE                                                                      \
void OpFunc                                                                   \
(                                                                             \
    FieldField<Field, ReturnType>& res,                                       \
    const FieldField<Field, Type1>& f                                         \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const FieldField<Field, Type1>& f                                         \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf                                  \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                       \
                                                                              \
TEMPLATE                                                                      \
void func                                                                     \
(                                                                             \
    FieldField<Field, ReturnType>& f,                                         \
    const FieldField<Field, Type1>& f1,                                       \
    const FieldField<Field, Type2>& f2                                        \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > func                                      \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const FieldField<Field, Type2>& f2                                        \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > func                                      \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > func                                      \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const FieldField<Field, Type2>& f2                                        \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > func                                      \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)               \
                                                                              \
TEMPLATE                                                                      \
void func                                                                     \
(                                                                             \
    FieldField<Field, ReturnType>& f,                                         \
    const Type1& s1,                                                          \
    const FieldField<Field, Type2>& f2                                        \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > func                                      \
(                                                                             \
    const Type1& s1,                                                          \
    const FieldField<Field, Type1>& f2                                        \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > func                                      \
(                                                                             \
    const Type1& s1,                                                          \
    const tmp<FieldField<Field, Type1> >& tf2                                 \
);


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)               \
                                                                              \
TEMPLATE                                                                      \
void func                                                                     \
(                                                                             \
    FieldField<Field, ReturnType>& f,                                         \
    const FieldField<Field, Type1>& f1,                                       \
    const Type2& s                                                            \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > func                                      \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const Type2& s                                                            \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > func                                      \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const Type2& s                                                            \
);


#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                  \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)                 \
                                                                              \
TEMPLATE                                                                      \
void OpFunc                                                                   \
(                                                                             \
    FieldField<Field, ReturnType>& f,                                         \
    const FieldField<Field, Type1>& f1,                                       \
    const FieldField<Field, Type2>& f2                                        \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const FieldField<Field, Type2>& f2                                        \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const FieldField<Field, Type2>& f2                                        \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)         \
                                                                              \
TEMPLATE                                                                      \
void OpFunc                                                                   \
(                                                                             \
    FieldField<Field, ReturnType>& f,                                         \
    const Type1& s1,                                                          \
    const FieldField<Field, Type2>& f2                                        \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const Type1& s1,                                                          \
    const FieldField<Field, Type2>& f2                                        \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const Type1& s1,                                                          \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
);


#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)         \
                                                                              \
TEMPLATE                                                                      \
void OpFunc                                                                   \
(                                                                             \
    FieldField<Field, ReturnType>& f,                                         \
    const FieldField<Field, Type1>& f1,                                       \
    const Type2& s2                                                           \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const Type2& s2                                                           \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<FieldField<Field, ReturnType> > operator Op                               \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const Type2& s2                                                           \
);


#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)            \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)


// ************************************************************************* //
