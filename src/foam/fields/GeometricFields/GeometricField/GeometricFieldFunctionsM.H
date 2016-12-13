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

#define UNARY_FUNCTION(ReturnType, Type1, Func, Dfunc)                        \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > Func                    \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& df                      \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > Func                    \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tdf1              \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_OPERATOR(ReturnType, Type1, Op, opFunc, Dfunc)                  \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > operator Op             \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& df1                     \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > operator Op             \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tdf1              \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                       \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > Func                    \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& df1,                    \
    const GeometricField<Type2, PatchField, GeoMesh>& df2                     \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > Func                    \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& df1,                    \
    const tmp<GeometricField<Type2, PatchField, GeoMesh> >& tdf2              \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > Func                    \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tdf1,             \
    const GeometricField<Type2, PatchField, GeoMesh>& df2                     \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > Func                    \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tdf1,             \
    const tmp<GeometricField<Type2, PatchField, GeoMesh> >& tdf2              \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)               \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > Func                    \
(                                                                             \
    const dimensioned<Type1>& dt1,                                            \
    const GeometricField<Type2, PatchField, GeoMesh>& df2                     \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > Func                    \
(                                                                             \
    const Type1& t1,                                                          \
    const GeometricField<Type2, PatchField, GeoMesh>& df2                     \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > Func                    \
(                                                                             \
    const dimensioned<Type1>& dt1,                                            \
    const tmp<GeometricField<Type2, PatchField, GeoMesh> >& tdf2              \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > Func                    \
(                                                                             \
    const Type1& t1,                                                          \
    const tmp<GeometricField<Type2, PatchField, GeoMesh> >& tdf2              \
);


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)               \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > Func                    \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& df1,                    \
    const dimensioned<Type2>& dt2                                             \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > Func                    \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& df1,                    \
    const Type2& t2                                                           \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > Func                    \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tdf1,             \
    const dimensioned<Type2>& dt2                                             \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > Func                    \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tdf2,             \
    const Type2& t2                                                           \
);


#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                  \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpName, OpFunc)         \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > operator Op             \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& df1,                    \
    const GeometricField<Type2, PatchField, GeoMesh>& df2                     \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > operator Op             \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& df1,                    \
    const tmp<GeometricField<Type2, PatchField, GeoMesh> >& tdf2              \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > operator Op             \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tdf1,             \
    const GeometricField<Type2, PatchField, GeoMesh>& df2                     \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > operator Op             \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tdf1,             \
    const tmp<GeometricField<Type2, PatchField, GeoMesh> >& tdf2              \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpName, OpFunc) \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > operator Op             \
(                                                                             \
    const dimensioned<Type1>& dt1,                                            \
    const GeometricField<Type2, PatchField, GeoMesh>& df2                     \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > operator Op             \
(                                                                             \
    const Type1& t1,                                                          \
    const GeometricField<Type2, PatchField, GeoMesh>& df2                     \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > operator Op             \
(                                                                             \
    const dimensioned<Type1>& dt1,                                            \
    const tmp<GeometricField<Type2, PatchField, GeoMesh> >& tdf2              \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > operator Op             \
(                                                                             \
    const Type1& t1,                                                          \
    const tmp<GeometricField<Type2, PatchField, GeoMesh> >& tdf2              \
);


#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpName, OpFunc) \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > operator Op             \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& df1,                    \
    const dimensioned<Type2>& dt2                                             \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > operator Op             \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& df1,                    \
    const Type2& t2                                                           \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > operator Op             \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tdf1,             \
    const dimensioned<Type2>& dt2                                             \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<ReturnType, PatchField, GeoMesh> > operator Op             \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tdf1,             \
    const Type2& t2                                                           \
);


#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpName, OpFunc)    \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpName, OpFunc)     \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpName, OpFunc)


// ************************************************************************* //
