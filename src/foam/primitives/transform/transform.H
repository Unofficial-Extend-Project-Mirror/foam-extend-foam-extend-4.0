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

InNamespace
    Foam

Description
    3D tensor transformation operations.

\*---------------------------------------------------------------------------*/

#ifndef transform_H
#define transform_H

#include "tensorField.H"
#include "symmTensor4thOrder.H"
#include "diagTensor.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline tensor rotationTensor
(
    const vector& n1,
    const vector& n2
)
{
    return
        (n1 & n2)*I
      + (1 - (n1 & n2))*sqr(n1 ^ n2)/(magSqr(n1 ^ n2) + VSMALL)
      + (n2*n1 - n1*n2);
}


inline tmp<tensorField> rotationTensor
(
    const vectorField& n1,
    const vectorField& n2
)
{
    return
        (n1 & n2)*I
      + (1 - (n1 & n2))*sqr(n1 ^ n2)/(magSqr(n1 ^ n2) + VSMALL)
      + (n2*n1 - n1*n2);
}


inline label transform(const tensor&, const label i)
{
    return i;
}


inline scalar transform(const tensor&, const scalar s)
{
    return s;
}


template<class Cmpt>
inline Vector<Cmpt> transform(const tensor& tt, const Vector<Cmpt>& v)
{
    return tt & v;
}


template<class Cmpt>
inline Tensor<Cmpt> transform(const tensor& tt, const Tensor<Cmpt>& t)
{
    //return tt & t & tt.T();
    return Tensor<Cmpt>
    (
        (tt.xx()*t.xx() + tt.xy()*t.yx() + tt.xz()*t.zx())*tt.xx()
      + (tt.xx()*t.xy() + tt.xy()*t.yy() + tt.xz()*t.zy())*tt.xy()
      + (tt.xx()*t.xz() + tt.xy()*t.yz() + tt.xz()*t.zz())*tt.xz(),

        (tt.xx()*t.xx() + tt.xy()*t.yx() + tt.xz()*t.zx())*tt.yx()
      + (tt.xx()*t.xy() + tt.xy()*t.yy() + tt.xz()*t.zy())*tt.yy()
      + (tt.xx()*t.xz() + tt.xy()*t.yz() + tt.xz()*t.zz())*tt.yz(),

        (tt.xx()*t.xx() + tt.xy()*t.yx() + tt.xz()*t.zx())*tt.zx()
      + (tt.xx()*t.xy() + tt.xy()*t.yy() + tt.xz()*t.zy())*tt.zy()
      + (tt.xx()*t.xz() + tt.xy()*t.yz() + tt.xz()*t.zz())*tt.zz(),

        (tt.yx()*t.xx() + tt.yy()*t.yx() + tt.yz()*t.zx())*tt.xx()
      + (tt.yx()*t.xy() + tt.yy()*t.yy() + tt.yz()*t.zy())*tt.xy()
      + (tt.yx()*t.xz() + tt.yy()*t.yz() + tt.yz()*t.zz())*tt.xz(),

        (tt.yx()*t.xx() + tt.yy()*t.yx() + tt.yz()*t.zx())*tt.yx()
      + (tt.yx()*t.xy() + tt.yy()*t.yy() + tt.yz()*t.zy())*tt.yy()
      + (tt.yx()*t.xz() + tt.yy()*t.yz() + tt.yz()*t.zz())*tt.yz(),

        (tt.yx()*t.xx() + tt.yy()*t.yx() + tt.yz()*t.zx())*tt.zx()
      + (tt.yx()*t.xy() + tt.yy()*t.yy() + tt.yz()*t.zy())*tt.zy()
      + (tt.yx()*t.xz() + tt.yy()*t.yz() + tt.yz()*t.zz())*tt.zz(),

        (tt.zx()*t.xx() + tt.zy()*t.yx() + tt.zz()*t.zx())*tt.xx()
      + (tt.zx()*t.xy() + tt.zy()*t.yy() + tt.zz()*t.zy())*tt.xy()
      + (tt.zx()*t.xz() + tt.zy()*t.yz() + tt.zz()*t.zz())*tt.xz(),

        (tt.zx()*t.xx() + tt.zy()*t.yx() + tt.zz()*t.zx())*tt.yx()
      + (tt.zx()*t.xy() + tt.zy()*t.yy() + tt.zz()*t.zy())*tt.yy()
      + (tt.zx()*t.xz() + tt.zy()*t.yz() + tt.zz()*t.zz())*tt.yz(),

        (tt.zx()*t.xx() + tt.zy()*t.yx() + tt.zz()*t.zx())*tt.zx()
      + (tt.zx()*t.xy() + tt.zy()*t.yy() + tt.zz()*t.zy())*tt.zy()
      + (tt.zx()*t.xz() + tt.zy()*t.yz() + tt.zz()*t.zz())*tt.zz()
    );
}


template<class Cmpt>
inline SphericalTensor<Cmpt> transform
(
    const tensor& tt,
    const SphericalTensor<Cmpt>& st
)
{
    return st;
}


template<class Cmpt>
inline SymmTensor4thOrder<Cmpt> transform
(
    const tensor& L,
    const SymmTensor4thOrder<Cmpt>& C
)
{
    //- represent fourth order tensors in 6x6 notation.  Rotation is given by
    //- C_rotated_af = A_ba * C_cd * A_ef
    //- where A is a function of L

    const scalar s = ::sqrt(2);
    const scalar A[6][6] =
    {
        {
            L.xx()*L.xx(), L.xy()*L.xy(), L.xz()*L.xz(), s*L.xx()*L.xy(), s*L.xy()*L.xz(), s*L.xz()*L.xx() },
        { L.yx()*L.yx(), L.yy()*L.yy(), L.yz()*L.yz(), s*L.yx()*L.yy(), s*L.yy()*L.yz(), s*L.yz()*L.yx() },
        { L.zx()*L.zx(), L.zy()*L.zy(), L.zz()*L.zz(), s*L.zx()*L.zy(), s*L.zy()*L.zz(), s*L.zz()*L.zx() },
        { s*L.xx()*L.yx(), s*L.xy()*L.yy(), s*L.xz()*L.yz(), (L.xy()*L.yx()+L.xx()*L.yy()), (L.xz()*L.yy()+L.xy()*L.yz()), (L.xx()*L.yz()+L.xz()*L.yx()) },
        { s*L.yx()*L.zx(), s*L.yy()*L.zy(), s*L.yz()*L.zz(), (L.yy()*L.zx()+L.yx()*L.zy()), (L.yz()*L.zy()+L.yy()*L.zz()), (L.yx()*L.zz()+L.yz()*L.zx()) },
        { s*L.zx()*L.xx(), s*L.zy()*L.xy(), s*L.zz()*L.xz(), (L.zy()*L.xx()+L.zx()*L.xy()), (L.zz()*L.xy()+L.zy()*L.xz()), (L.zx()*L.xz()+L.zz()*L.xx()) }
    };

    return symmTensor4thOrder
    (
        // xxxx
        A[0][0] * C.xxxx() * A[0][0] +
        A[1][0] * C.xxyy() * A[0][0] +
        A[2][0] * C.xxzz() * A[0][0] +
        A[0][0] * C.xxyy() * A[1][0] +
        A[1][0] * C.yyyy() * A[1][0] +
        A[2][0] * C.yyzz() * A[1][0] +
        A[0][0] * C.xxzz() * A[2][0] +
        A[1][0] * C.yyzz() * A[2][0] +
        A[2][0] * C.zzzz() * A[2][0] +
        A[3][0] * C.xyxy() * A[3][0] +
        A[4][0] * C.yzyz() * A[4][0] +
        A[5][0] * C.zxzx() * A[5][0],

        // xxyy
        A[0][0] * C.xxxx() * A[0][1] +
        A[1][0] * C.xxyy() * A[0][1] +
        A[2][0] * C.xxzz() * A[0][1] +
        A[0][0] * C.xxyy() * A[1][1] +
        A[1][0] * C.yyyy() * A[1][1] +
        A[2][0] * C.yyzz() * A[1][1] +
        A[0][0] * C.xxzz() * A[2][1] +
        A[1][0] * C.yyzz() * A[2][1] +
        A[2][0] * C.zzzz() * A[2][1] +
        A[3][0] * C.xyxy() * A[3][1] +
        A[4][0] * C.yzyz() * A[4][1] +
        A[5][0] * C.zxzx() * A[5][1],

        // xzz
        A[0][0] * C.xxxx() * A[0][2] +
        A[1][0] * C.xxyy() * A[0][2] +
        A[2][0] * C.xxzz() * A[0][2] +
        A[0][0] * C.xxyy() * A[1][2] +
        A[1][0] * C.yyyy() * A[1][2] +
        A[2][0] * C.yyzz() * A[1][2] +
        A[0][0] * C.xxzz() * A[2][2] +
        A[1][0] * C.yyzz() * A[2][2] +
        A[2][0] * C.zzzz() * A[2][2] +
        A[3][0] * C.xyxy() * A[3][2] +
        A[4][0] * C.yzyz() * A[4][2] +
        A[5][0] * C.zxzx() * A[5][2],

        // yyyy
        A[0][1] * C.xxxx() * A[0][1] +
        A[1][1] * C.xxyy() * A[0][1] +
        A[2][1] * C.xxzz() * A[0][1] +
        A[0][1] * C.xxyy() * A[1][1] +
        A[1][1] * C.yyyy() * A[1][1] +
        A[2][1] * C.yyzz() * A[1][1] +
        A[0][1] * C.xxzz() * A[2][1] +
        A[1][1] * C.yyzz() * A[2][1] +
        A[2][1] * C.zzzz() * A[2][1] +
        A[3][1] * C.xyxy() * A[3][1] +
        A[4][1] * C.yzyz() * A[4][1] +
        A[5][1] * C.zxzx() * A[5][1],

        // yyzz
        A[0][1] * C.xxxx() * A[0][2] +
        A[1][1] * C.xxyy() * A[0][2] +
        A[2][1] * C.xxzz() * A[0][2] +
        A[0][1] * C.xxyy() * A[1][2] +
        A[1][1] * C.yyyy() * A[1][2] +
        A[2][1] * C.yyzz() * A[1][2] +
        A[0][1] * C.xxzz() * A[2][2] +
        A[1][1] * C.yyzz() * A[2][2] +
        A[2][1] * C.zzzz() * A[2][2] +
        A[3][1] * C.xyxy() * A[3][2] +
        A[4][1] * C.yzyz() * A[4][2] +
        A[5][1] * C.zxzx() * A[5][2],

        // zzzz
        A[0][2] * C.xxxx() * A[0][2] +
        A[1][2] * C.xxyy() * A[0][2] +
        A[2][2] * C.xxzz() * A[0][2] +
        A[0][2] * C.xxyy() * A[1][2] +
        A[1][2] * C.yyyy() * A[1][2] +
        A[2][2] * C.yyzz() * A[1][2] +
        A[0][2] * C.xxzz() * A[2][2] +
        A[1][2] * C.yyzz() * A[2][2] +
        A[2][2] * C.zzzz() * A[2][2] +
        A[3][2] * C.xyxy() * A[3][2] +
        A[4][2] * C.yzyz() * A[4][2] +
        A[5][2] * C.zxzx() * A[5][2],

        // xyxy
        A[0][3] * C.xxxx() * A[0][3] +
        A[1][3] * C.xxyy() * A[0][3] +
        A[2][3] * C.xxzz() * A[0][3] +
        A[0][3] * C.xxyy() * A[1][3] +
        A[1][3] * C.yyyy() * A[1][3] +
        A[2][3] * C.yyzz() * A[1][3] +
        A[0][3] * C.xxzz() * A[2][3] +
        A[1][3] * C.yyzz() * A[2][3] +
        A[2][3] * C.zzzz() * A[2][3] +
        A[3][3] * C.xyxy() * A[3][3] +
        A[4][3] * C.yzyz() * A[4][3] +
        A[5][3] * C.zxzx() * A[5][3],

        // yzyz
        A[0][4] * C.xxxx() * A[0][4] +
        A[1][4] * C.xxyy() * A[0][4] +
        A[2][4] * C.xxzz() * A[0][4] +
        A[0][4] * C.xxyy() * A[1][4] +
        A[1][4] * C.yyyy() * A[1][4] +
        A[2][4] * C.yyzz() * A[1][4] +
        A[0][4] * C.xxzz() * A[2][4] +
        A[1][4] * C.yyzz() * A[2][4] +
        A[2][4] * C.zzzz() * A[2][4] +
        A[3][4] * C.xyxy() * A[3][4] +
        A[4][4] * C.yzyz() * A[4][4] +
        A[5][4] * C.zxzx() * A[5][4],

        // zxzx
        A[0][5] * C.xxxx() * A[0][5] +
        A[1][5] * C.xxyy() * A[0][5] +
        A[2][5] * C.xxzz() * A[0][5] +
        A[0][5] * C.xxyy() * A[1][5] +
        A[1][5] * C.yyyy() * A[1][5] +
        A[2][5] * C.yyzz() * A[1][5] +
        A[0][5] * C.xxzz() * A[2][5] +
        A[1][5] * C.yyzz() * A[2][5] +
        A[2][5] * C.zzzz() * A[2][5] +
        A[3][5] * C.xyxy() * A[3][5] +
        A[4][5] * C.yzyz() * A[4][5] +
        A[5][5] * C.zxzx() * A[5][5]
    );
}


template<class Cmpt>
inline DiagTensor<Cmpt> transform
(
    const tensor& tt,
    const DiagTensor<Cmpt>& st
)
{
    notImplemented
    (
        "transform.H\n"
        "template<>\n"
        "inline DiagTensor<Cmpt> transform\n"
        "(\n"
        "const tensor& tt,\n"
        "const DiagTensor<Cmpt>& st\n"
        ")\n"
        "not implemented"
    );

    return st;
}


template<class Cmpt>
inline SymmTensor<Cmpt> transform(const tensor& tt, const SymmTensor<Cmpt>& st)
{
    return SymmTensor<Cmpt>
    (
        (tt.xx()*st.xx() + tt.xy()*st.xy() + tt.xz()*st.xz())*tt.xx()
      + (tt.xx()*st.xy() + tt.xy()*st.yy() + tt.xz()*st.yz())*tt.xy()
      + (tt.xx()*st.xz() + tt.xy()*st.yz() + tt.xz()*st.zz())*tt.xz(),

        (tt.xx()*st.xx() + tt.xy()*st.xy() + tt.xz()*st.xz())*tt.yx()
      + (tt.xx()*st.xy() + tt.xy()*st.yy() + tt.xz()*st.yz())*tt.yy()
      + (tt.xx()*st.xz() + tt.xy()*st.yz() + tt.xz()*st.zz())*tt.yz(),

        (tt.xx()*st.xx() + tt.xy()*st.xy() + tt.xz()*st.xz())*tt.zx()
      + (tt.xx()*st.xy() + tt.xy()*st.yy() + tt.xz()*st.yz())*tt.zy()
      + (tt.xx()*st.xz() + tt.xy()*st.yz() + tt.xz()*st.zz())*tt.zz(),

        (tt.yx()*st.xx() + tt.yy()*st.xy() + tt.yz()*st.xz())*tt.yx()
      + (tt.yx()*st.xy() + tt.yy()*st.yy() + tt.yz()*st.yz())*tt.yy()
      + (tt.yx()*st.xz() + tt.yy()*st.yz() + tt.yz()*st.zz())*tt.yz(),

        (tt.yx()*st.xx() + tt.yy()*st.xy() + tt.yz()*st.xz())*tt.zx()
      + (tt.yx()*st.xy() + tt.yy()*st.yy() + tt.yz()*st.yz())*tt.zy()
      + (tt.yx()*st.xz() + tt.yy()*st.yz() + tt.yz()*st.zz())*tt.zz(),

        (tt.zx()*st.xx() + tt.zy()*st.xy() + tt.zz()*st.xz())*tt.zx()
      + (tt.zx()*st.xy() + tt.zy()*st.yy() + tt.zz()*st.yz())*tt.zy()
      + (tt.zx()*st.xz() + tt.zy()*st.yz() + tt.zz()*st.zz())*tt.zz()
    );
}


template<class Type1, class Type2>
inline Type1 transformMask(const Type2& t)
{
    return t;
}


template<>
inline sphericalTensor transformMask<sphericalTensor>(const tensor& t)
{
    return sph(t);
}


template<>
inline symmTensor transformMask<symmTensor>(const tensor& t)
{
    return symm(t);
}


template<>
inline symmTensor4thOrder transformMask<symmTensor4thOrder>(const tensor& t)
{
  notImplemented("transform.H\n"
         "template<>\n"
         "inline symmTensor4thOrder transformMask<symmTensor4thOrder>(const tensor& t)\n"
         "not implemented");

  return symmTensor4thOrder::zero;
}


template<>
inline diagTensor transformMask<diagTensor>(const tensor& t)
{
  return diagTensor(t.xx(), t.yy(), t.zz());
}


//- Estimate angle of vec in coordinate system (e0, e1, e0^e1).
//  Is guaranteed to return increasing number but is not correct
//  angle. Used for sorting angles.  All input vectors need to be normalized.
//
// Calculates scalar which increases with angle going from e0 to vec in
// the coordinate system e0, e1, e0^e1
//
// Jumps from 2PI -> 0 at -SMALL so parallel vectors with small rounding errors
// should hopefully still get the same quadrant.
//
inline scalar pseudoAngle
(
    const vector& e0,
    const vector& e1,
    const vector& vec
)
{
    scalar cos = vec & e0;
    scalar sin = vec & e1;

    if (sin < -SMALL)
    {
        return (3.0 + cos)*mathematicalConstant::piByTwo;
    }
    else
    {
        return (1.0 - cos)*mathematicalConstant::piByTwo;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
